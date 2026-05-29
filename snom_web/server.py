from __future__ import annotations

import argparse
import json
import mimetypes
import multiprocessing
import queue as queue_module
import threading
import time
import uuid
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

from .calculator import ValidationError, calculate_spectrum
from .materials import list_materials


PACKAGE_ROOT = Path(__file__).resolve().parent
STATIC_ROOT = PACKAGE_ROOT / "static"
TIP_SCHEME_PATHS = {
    "bulk": PACKAGE_ROOT.parent / "Tip_schemes" / "Tip_scheme.jpg",
    "layered": PACKAGE_ROOT.parent / "Tip_schemes" / "Tip_scheme_layers.jpg",
}
TIMEOUT_SECONDS = 600

ACTIVE_CALCULATIONS: dict[str, multiprocessing.Process] = {}
CANCELLED_CALCULATIONS: set[str] = set()
ACTIVE_CALCULATIONS_LOCK = threading.Lock()

EXAMPLE_PERMITTIVITY_FILE = """# frequency_cm-1 Re_epsilon Im_epsilon
1650 2.22 0.02
1675 2.24 0.03
1700 2.28 0.04
1725 2.32 0.06
1750 2.35 0.08
1775 2.33 0.07
1800 2.30 0.05
"""

EXAMPLE_EXPERIMENTAL_FILE = """# frequency_cm-1 amplitude phase
1650 0.93 -0.55
1675 0.98 -0.51
1700 1.04 -0.46
1725 1.08 -0.40
1750 1.12 -0.33
1775 1.10 -0.29
1800 1.06 -0.25
"""


class AppHandler(BaseHTTPRequestHandler):
    server_version = "sSNOMWeb/1.0"

    def do_GET(self) -> None:
        if self.path == "/api/materials":
            self._send_json(HTTPStatus.OK, {"materials": list_materials()})
            return

        if self.path == "/api/example/permittivity":
            self._send_download("permittivity-example.txt", EXAMPLE_PERMITTIVITY_FILE)
            return

        if self.path == "/api/example/experimental":
            self._send_download("experimental-spectrum-example.txt", EXAMPLE_EXPERIMENTAL_FILE)
            return

        if self.path in {"/assets/tip-scheme.jpg", "/assets/tip-scheme-bulk.jpg"}:
            self._send_file(TIP_SCHEME_PATHS["bulk"])
            return

        if self.path == "/assets/tip-scheme-layered.jpg":
            self._send_file(TIP_SCHEME_PATHS["layered"])
            return

        if self.path in {"/", "/index.html"}:
            self._send_file(STATIC_ROOT / "index.html")
            return

        if self.path.startswith("/static/"):
            relative_path = self.path.removeprefix("/static/")
            static_path = (STATIC_ROOT / relative_path).resolve()
            if STATIC_ROOT.resolve() not in static_path.parents and static_path != STATIC_ROOT.resolve():
                self._send_json(HTTPStatus.NOT_FOUND, {"error": "File not found."})
                return
            self._send_file(static_path)
            return

        self._send_json(HTTPStatus.NOT_FOUND, {"error": "Route not found."})

    def do_POST(self) -> None:
        if self.path == "/api/cancel":
            self._handle_cancel()
            return

        if self.path != "/api/calculate":
            self._send_json(HTTPStatus.NOT_FOUND, {"error": "Route not found."})
            return

        try:
            payload = self._read_json_body()
        except json.JSONDecodeError:
            self._send_json(HTTPStatus.BAD_REQUEST, {"error": "Request body must contain valid JSON."})
            return

        calculation_id = str(payload.get("calculationId") or uuid.uuid4())
        try:
            result = _run_calculation_with_timeout(payload, calculation_id)
        except ValidationError as error:
            self._send_json(HTTPStatus.BAD_REQUEST, {"error": str(error)})
            return
        except CalculationCancelled as error:
            self._send_json(HTTPStatus.CONFLICT, {"error": str(error)})
            return
        except TimeoutError as error:
            self._send_json(HTTPStatus.GATEWAY_TIMEOUT, {"error": str(error)})
            return
        except Exception as error:  # pragma: no cover - protects the HTTP layer
            self._send_json(HTTPStatus.INTERNAL_SERVER_ERROR, {"error": str(error)})
            return

        self._send_json(HTTPStatus.OK, result)

    def _handle_cancel(self) -> None:
        try:
            payload = self._read_json_body()
        except json.JSONDecodeError:
            self._send_json(HTTPStatus.BAD_REQUEST, {"error": "Request body must contain valid JSON."})
            return

        calculation_id = payload.get("calculationId")
        if not isinstance(calculation_id, str) or not calculation_id:
            self._send_json(HTTPStatus.BAD_REQUEST, {"error": "calculationId is required."})
            return

        if cancel_calculation(calculation_id):
            self._send_json(HTTPStatus.OK, {"status": "cancelled"})
            return

        self._send_json(HTTPStatus.NOT_FOUND, {"error": "Calculation not found."})

    def _read_json_body(self) -> dict[str, object]:
        content_length = int(self.headers.get("Content-Length", "0"))
        payload = json.loads(self.rfile.read(content_length).decode("utf-8"))
        if not isinstance(payload, dict):
            raise json.JSONDecodeError("JSON body must be an object.", "", 0)
        return payload

    def _send_json(self, status: HTTPStatus, payload: dict[str, object]) -> None:
        body = json.dumps(payload, ensure_ascii=False).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def _send_download(self, filename: str, content: str) -> None:
        body = content.encode("utf-8")
        self.send_response(HTTPStatus.OK)
        self.send_header("Content-Type", "text/plain; charset=utf-8")
        self.send_header("Content-Disposition", f'attachment; filename="{filename}"')
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def _send_file(self, path: Path) -> None:
        if not path.exists() or not path.is_file():
            self._send_json(HTTPStatus.NOT_FOUND, {"error": "File not found."})
            return

        content_type, _ = mimetypes.guess_type(path.name)
        body = path.read_bytes()
        self.send_response(HTTPStatus.OK)
        self.send_header("Content-Type", content_type or "application/octet-stream")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)


class CalculationCancelled(RuntimeError):
    """Raised when a user explicitly stops an active calculation."""


def cancel_calculation(calculation_id: str) -> bool:
    with ACTIVE_CALCULATIONS_LOCK:
        process = ACTIVE_CALCULATIONS.get(calculation_id)
        if process is None:
            return False
        CANCELLED_CALCULATIONS.add(calculation_id)

    if process.is_alive():
        process.terminate()
        process.join(timeout=1.0)
        if process.is_alive():
            process.kill()
            process.join(timeout=1.0)

    return True


def _run_calculation_with_timeout(
    payload: dict[str, object], calculation_id: str
) -> dict[str, object]:
    context = multiprocessing.get_context("spawn")
    result_queue: multiprocessing.Queue[dict[str, object]] = context.Queue()
    process = context.Process(target=_calculation_worker, args=(payload, result_queue))
    process.start()
    deadline = time.monotonic() + TIMEOUT_SECONDS

    with ACTIVE_CALCULATIONS_LOCK:
        ACTIVE_CALCULATIONS[calculation_id] = process

    try:
        while True:
            if time.monotonic() >= deadline:
                cancel_calculation(calculation_id)
                raise TimeoutError("Time out.")

            try:
                message = result_queue.get(timeout=0.1)
            except queue_module.Empty:
                if process.is_alive():
                    continue
                if _is_calculation_cancelled(calculation_id):
                    raise CalculationCancelled("Calculation stopped.")
                try:
                    message = result_queue.get_nowait()
                except queue_module.Empty:
                    if process.exitcode == 0:
                        raise RuntimeError("The calculation finished without returning a result.")
                    raise RuntimeError("The calculation process stopped unexpectedly.")

            if _is_calculation_cancelled(calculation_id):
                raise CalculationCancelled("Calculation stopped.")

            process.join(timeout=1.0)
            if process.is_alive():
                process.terminate()
                process.join(timeout=1.0)
            return _resolve_worker_message(message)
    finally:
        with ACTIVE_CALCULATIONS_LOCK:
            ACTIVE_CALCULATIONS.pop(calculation_id, None)
            CANCELLED_CALCULATIONS.discard(calculation_id)
        result_queue.close()
        result_queue.join_thread()


def _is_calculation_cancelled(calculation_id: str) -> bool:
    with ACTIVE_CALCULATIONS_LOCK:
        return calculation_id in CANCELLED_CALCULATIONS


def _resolve_worker_message(message: dict[str, object]) -> dict[str, object]:
    if message["status"] == "validation_error":
        raise ValidationError(str(message["error"]))
    if message["status"] == "error":
        raise RuntimeError(str(message["error"]))
    return message["result"]


def _calculation_worker(
    payload: dict[str, object], queue: multiprocessing.Queue[dict[str, object]]
) -> None:
    try:
        queue.put({"status": "ok", "result": calculate_spectrum(payload)})
    except ValidationError as error:
        queue.put({"status": "validation_error", "error": str(error)})
    except Exception as error:  # pragma: no cover - child process guard
        queue.put({"status": "error", "error": str(error)})


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the s-SNOM web calculator.")
    parser.add_argument("--host", default="127.0.0.1", help="Host to bind to.")
    parser.add_argument("--port", type=int, default=8000, help="Port to bind to.")
    arguments = parser.parse_args()

    server = ThreadingHTTPServer((arguments.host, arguments.port), AppHandler)
    print(f"s-SNOM web calculator is running at http://{arguments.host}:{arguments.port}")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()
