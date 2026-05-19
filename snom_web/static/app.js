const state = {
  latestResult: null,
};

const form = document.getElementById("calculator-form");
const calculateButton = document.getElementById("calculate-button");
const downloadResultsButton = document.getElementById("download-results-button");
const formStatus = document.getElementById("form-status");
const errorBox = document.getElementById("error-box");
const builtInMaterialSelect = document.getElementById("built-in-material");
const sampleBuiltIn = document.getElementById("sample-built-in");
const sampleUpload = document.getElementById("sample-upload");
const farFieldEnabled = document.getElementById("far-field-enabled");
const farFieldFields = document.getElementById("far-field-fields");
const experimentalEnabled = document.getElementById("experimental-enabled");
const experimentalFields = document.getElementById("experimental-fields");
const experimentalFile = document.getElementById("experimental-file");
const discretizationFields = document.getElementById("discretization-fields");
const resultSummary = document.getElementById("result-summary");
const chartModal = document.getElementById("chart-modal");
const chartModalTitle = document.getElementById("chart-modal-title");
const chartModalPlot = document.getElementById("chart-modal-plot");
const chartModalClose = document.getElementById("chart-modal-close");


async function initialize() {
  await loadMaterials();
  syncUiState();

  form.addEventListener("change", syncUiState);
  form.addEventListener("submit", handleSubmit);
  downloadResultsButton.addEventListener("click", handleDownloadResults);
  chartModalClose.addEventListener("click", closeChartModal);
  chartModal.addEventListener("click", (event) => {
    if (event.target === chartModal) {
      closeChartModal();
    }
  });
  document.addEventListener("keydown", (event) => {
    if (event.key === "Escape" && !chartModal.classList.contains("hidden")) {
      closeChartModal();
    }
  });
  document.querySelectorAll(".plot").forEach((plot) => {
    plot.addEventListener("click", () => openChartModal(plot.id));
  });
  document.querySelectorAll(".chart-open-button").forEach((button) => {
    button.addEventListener("click", () => openChartModal(button.dataset.chartTarget));
  });

  renderEmptyPlots();
}


async function loadMaterials() {
  const response = await fetch("/api/materials");
  const payload = await response.json();

  builtInMaterialSelect.innerHTML = "";
  for (const material of payload.materials) {
    const option = document.createElement("option");
    option.value = material.id;
    option.textContent = material.label;
    builtInMaterialSelect.appendChild(option);
  }
}


function syncUiState() {
  const sampleInputMethod = getCheckedValue("sampleInputMethod");
  const algorithm = getCheckedValue("algorithm");

  document.querySelectorAll(".choice-card").forEach((card) => {
    const input = card.querySelector("input[type='radio']");
    card.classList.toggle("active", Boolean(input && input.checked));
  });

  sampleBuiltIn.classList.toggle("hidden", sampleInputMethod !== "builtIn");
  sampleUpload.classList.toggle("hidden", sampleInputMethod !== "upload");

  farFieldFields.classList.toggle("is-disabled", !farFieldEnabled.checked);
  farFieldFields.querySelectorAll("input").forEach((input) => {
    input.disabled = !farFieldEnabled.checked;
  });

  experimentalFields.classList.toggle("is-disabled", !experimentalEnabled.checked);
  experimentalFile.disabled = !experimentalEnabled.checked;

  const accurateSelected = algorithm === "accurate";
  discretizationFields.classList.toggle("is-disabled", !accurateSelected);
  document.getElementById("discretization-n").disabled = !accurateSelected;
  document.getElementById("discretization-m").disabled = !accurateSelected;
}


async function handleSubmit(event) {
  event.preventDefault();
  clearError();
  setStatus("Preparing calculation request...");

  try {
    const payload = await buildPayload();
    const validationErrors = validatePayload(payload);
    if (validationErrors.length > 0) {
      showError(validationErrors.join("\n"));
      setStatus("Fix the validation errors and try again.");
      return;
    }

    setBusy(true, "Running calculation...");
    const response = await fetch("/api/calculate", {
      method: "POST",
      headers: {
        "Content-Type": "application/json",
      },
      body: JSON.stringify(payload),
    });
    const result = await response.json();

    if (!response.ok) {
      throw new Error(result.error || "Calculation failed.");
    }

    state.latestResult = result;
    renderResult(result);
    downloadResultsButton.hidden = false;
    setStatus("Calculation finished.");
  } catch (error) {
    showError(error.message || "Calculation failed.");
    setStatus("Calculation failed.");
  } finally {
    setBusy(false, "");
  }
}


async function buildPayload() {
  const sampleInputMethod = getCheckedValue("sampleInputMethod");
  const experimentalOn = experimentalEnabled.checked;

  return {
    modelType: getCheckedValue("modelType"),
    algorithm: getCheckedValue("algorithm"),
    tip: {
      radius: getNumber("tip-radius"),
      length: getNumber("tip-length"),
      amplitude: getNumber("tip-amplitude"),
      sampleGap: getNumber("tip-sample-gap"),
      referenceGap: getNumber("tip-reference-gap"),
      epsilonReal: getNumber("tip-epsilon-real"),
      epsilonImag: getNumber("tip-epsilon-imag"),
    },
    farField: {
      enabled: farFieldEnabled.checked,
      coefficient: getNumber("far-field-coefficient"),
      angleDegrees: getNumber("far-field-angle"),
    },
    reference: {
      epsilonReal: getNumber("reference-epsilon-real"),
      epsilonImag: getNumber("reference-epsilon-imag"),
    },
    harmonicNumber: getInteger("harmonic-number"),
    frequencyRange: {
      start: getNumber("frequency-start"),
      end: getNumber("frequency-end"),
      points: getInteger("frequency-points"),
    },
    sample: {
      inputMethod: sampleInputMethod,
      builtInMaterial: sampleInputMethod === "builtIn" ? builtInMaterialSelect.value : null,
      uploadedPermittivity:
        sampleInputMethod === "upload"
          ? await readFileContent(document.getElementById("permittivity-file"))
          : null,
    },
    experimentalSpectrum: {
      enabled: experimentalOn,
      content: experimentalOn ? await readFileContent(experimentalFile) : null,
    },
    discretization: {
      N: getInteger("discretization-n"),
      M: getInteger("discretization-m"),
    },
  };
}


function validatePayload(payload) {
  const errors = [];

  if (!(payload.tip.radius > 0)) {
    errors.push("R must be greater than 0.");
  }
  if (!(payload.tip.length > 0)) {
    errors.push("L must be greater than 0.");
  }
  if (!(payload.tip.length > 2 * payload.tip.radius)) {
    errors.push("L must be greater than 2R so the spheroid geometry remains valid.");
  }
  if (!(payload.tip.amplitude > 0)) {
    errors.push("A must be greater than 0.");
  }
  if (!(payload.tip.sampleGap > 0)) {
    errors.push("H0,s must be greater than 0.");
  }
  if (!(payload.tip.referenceGap > 0)) {
    errors.push("H0,r must be greater than 0.");
  }
  if (!(payload.harmonicNumber > 0)) {
    errors.push("Harmonic number must be a positive integer.");
  }
  if (!(payload.frequencyRange.start < payload.frequencyRange.end)) {
    errors.push("Start frequency must be smaller than end frequency.");
  }
  if (!(payload.frequencyRange.points >= 1 && payload.frequencyRange.points <= 1000)) {
    errors.push("Number of frequency points must be between 1 and 1000.");
  }
  if (payload.farField.angleDegrees < 0 || payload.farField.angleDegrees > 90) {
    errors.push("Illumination angle must be between 0 and 90 degrees.");
  }

  if (payload.algorithm === "accurate") {
    if (!(payload.discretization.N >= 1 && payload.discretization.N <= 500)) {
      errors.push("N must be between 1 and 500.");
    }
    if (!(payload.discretization.M >= 1 && payload.discretization.M <= 100)) {
      errors.push("M must be between 1 and 100.");
    }
  }

  if (payload.sample.inputMethod === "upload" && !payload.sample.uploadedPermittivity) {
    errors.push("Please upload a permittivity file.");
  }
  if (payload.experimentalSpectrum.enabled && !payload.experimentalSpectrum.content) {
    errors.push("Please upload an experimental spectrum file.");
  }

  return errors;
}


function renderResult(result) {
  updateSummary(result.metadata);
  renderPlots(result);
}


function updateSummary(metadata) {
  resultSummary.classList.remove("muted-summary");
  resultSummary.innerHTML = `
    <div>
      <span class="summary-label">State</span>
      <strong>Ready</strong>
    </div>
    <div>
      <span class="summary-label">Material</span>
      <strong>${escapeHtml(metadata.materialLabel)}</strong>
    </div>
    <div>
      <span class="summary-label">Algorithm</span>
      <strong>${metadata.algorithm === "fast" ? "Fast" : "Accurate"}</strong>
    </div>
    <div>
      <span class="summary-label">Harmonic</span>
      <strong>n = ${metadata.harmonicNumber}</strong>
    </div>
    <div>
      <span class="summary-label">Frequency grid</span>
      <strong>${metadata.frequencyStart} to ${metadata.frequencyEnd} cm-1</strong>
    </div>
    <div>
      <span class="summary-label">Points</span>
      <strong>${metadata.frequencyPoints}</strong>
    </div>
    <div>
      <span class="summary-label">Sample discretization</span>
      <strong>N = ${metadata.sampleDiscretizationN}, M = ${metadata.sampleDiscretizationM}</strong>
    </div>
    <div>
      <span class="summary-label">Far-field</span>
      <strong>${metadata.farFieldEnabled ? "Enabled" : "Disabled"}</strong>
    </div>
  `;
}


function renderPlots(result) {
  const x = result.frequency;
  const plotConfig = {
    responsive: true,
    displaylogo: false,
  };

  Plotly.react(
    "permittivity-plot",
    [
      {
        x,
        y: result.epsilonReal,
        name: "Re(ε)",
        mode: "lines",
        line: { color: "#3157d8", width: 2.5 },
      },
      {
        x,
        y: result.epsilonImag,
        name: "Im(ε)",
        mode: "lines",
        line: { color: "#d43f6c", width: 2.5 },
      },
    ],
    baseLayout("Dielectric permittivity, ε", "ω, cm-1", "ε"),
    plotConfig,
  );

  const amplitudeTraces = [
    {
      x,
      y: result.amplitude,
      name: "Calculated spectrum",
      mode: "lines",
      line: { color: "#3157d8", width: 2.5 },
    },
  ];
  const phaseTraces = [
    {
      x,
      y: result.phase,
      name: "Calculated spectrum",
      mode: "lines",
      line: { color: "#3157d8", width: 2.5 },
    },
  ];

  if (result.experimental) {
    amplitudeTraces.push({
      x: result.experimental.frequency,
      y: result.experimental.amplitude,
      name: "Experimental spectrum",
      mode: "lines+markers",
      marker: { size: 5, color: "#ff8f3f" },
      line: { color: "#ff8f3f", width: 1.8 },
    });
    phaseTraces.push({
      x: result.experimental.frequency,
      y: result.experimental.phase,
      name: "Experimental spectrum",
      mode: "lines+markers",
      marker: { size: 5, color: "#ff8f3f" },
      line: { color: "#ff8f3f", width: 1.8 },
    });
  }

  Plotly.react(
    "amplitude-plot",
    amplitudeTraces,
    baseLayout("Near-field amplitude, |σn|", "ω, cm-1", "|σn|"),
    plotConfig,
  );

  Plotly.react(
    "phase-plot",
    phaseTraces,
    baseLayout("Near-field phase, arg(σn)", "ω, cm-1", "arg(σn), rad"),
    plotConfig,
  );
}


function renderEmptyPlots() {
  const emptyLayout = baseLayout("Awaiting calculation", "ω, cm-1", "");
  emptyLayout.annotations = [
    {
      text: "Run a calculation to display the spectrum.",
      showarrow: false,
      xref: "paper",
      yref: "paper",
      x: 0.5,
      y: 0.5,
      font: { size: 15, color: "#5e6a83" },
    },
  ];
  Plotly.newPlot("permittivity-plot", [], emptyLayout, { responsive: true, displaylogo: false });
  Plotly.newPlot("amplitude-plot", [], emptyLayout, { responsive: true, displaylogo: false });
  Plotly.newPlot("phase-plot", [], emptyLayout, { responsive: true, displaylogo: false });
}


function baseLayout(title, xLabel, yLabel) {
  return {
    title: { text: title, font: { size: 12 } },
    paper_bgcolor: "rgba(0,0,0,0)",
    plot_bgcolor: "#ffffff",
    autosize: true,
    margin: { l: 42, r: 10, t: 30, b: 34 },
    xaxis: {
      title: xLabel,
      titlefont: { size: 10 },
      tickfont: { size: 9 },
      gridcolor: "#e2e8f4",
      zerolinecolor: "#e2e8f4",
    },
    yaxis: {
      title: yLabel,
      titlefont: { size: 10 },
      tickfont: { size: 9 },
      gridcolor: "#e2e8f4",
      zerolinecolor: "#e2e8f4",
    },
    legend: {
      orientation: "h",
      x: 0,
      y: 1.18,
      font: { size: 9 },
    },
  };
}


function openChartModal(chartId) {
  const source = document.getElementById(chartId);
  if (!source || !source.data || !source.layout) {
    return;
  }

  chartModalTitle.textContent = source.dataset.chartTitle || "Chart";
  chartModal.classList.remove("hidden");

  const data = JSON.parse(JSON.stringify(source.data));
  const layout = JSON.parse(JSON.stringify(source.layout));
  layout.autosize = true;
  layout.margin = { l: 78, r: 28, t: 62, b: 68 };
  layout.title = {
    text: chartModalTitle.textContent,
    font: { size: 20 },
  };
  layout.legend = {
    ...layout.legend,
    orientation: "h",
    x: 0,
    y: 1.12,
    font: { size: 13 },
  };
  layout.xaxis = {
    ...layout.xaxis,
    titlefont: { size: 14 },
    tickfont: { size: 12 },
  };
  layout.yaxis = {
    ...layout.yaxis,
    titlefont: { size: 14 },
    tickfont: { size: 12 },
  };

  Plotly.react(chartModalPlot, data, layout, { responsive: true, displaylogo: false });
  requestAnimationFrame(() => Plotly.Plots.resize(chartModalPlot));
}


function closeChartModal() {
  chartModal.classList.add("hidden");
  Plotly.purge(chartModalPlot);
}


function handleDownloadResults() {
  if (!state.latestResult) {
    return;
  }

  const rows = [
    [
      "frequency_cm-1",
      "Re_epsilon_sample",
      "Im_epsilon_sample",
      "Re_sigma_n",
      "Im_sigma_n",
      "abs_sigma_n",
      "arg_sigma_n",
    ],
  ];

  const result = state.latestResult;
  for (let index = 0; index < result.frequency.length; index += 1) {
    rows.push([
      result.frequency[index],
      result.epsilonReal[index],
      result.epsilonImag[index],
      result.sigmaReal[index],
      result.sigmaImag[index],
      result.amplitude[index],
      result.phase[index],
    ]);
  }

  const csvContent = rows.map((row) => row.join(",")).join("\n");
  const blob = new Blob([csvContent], { type: "text/csv;charset=utf-8" });
  const url = URL.createObjectURL(blob);
  const anchor = document.createElement("a");
  anchor.href = url;
  anchor.download = result.suggestedFilename || "s-snom-results.csv";
  document.body.appendChild(anchor);
  anchor.click();
  anchor.remove();
  URL.revokeObjectURL(url);
}


function setBusy(isBusy, statusText) {
  calculateButton.disabled = isBusy;
  calculateButton.textContent = isBusy ? "Calculating..." : "Calculate";
  if (statusText) {
    setStatus(statusText);
  }
}


function setStatus(message) {
  formStatus.textContent = message;
}


function showError(message) {
  errorBox.textContent = message;
  errorBox.classList.remove("hidden");
}


function clearError() {
  errorBox.textContent = "";
  errorBox.classList.add("hidden");
}


async function readFileContent(input) {
  const file = input.files && input.files[0];
  if (!file) {
    return null;
  }
  return file.text();
}


function getNumber(id) {
  return Number(document.getElementById(id).value);
}


function getInteger(id) {
  return Number.parseInt(document.getElementById(id).value, 10);
}


function getCheckedValue(name) {
  return document.querySelector(`input[name="${name}"]:checked`).value;
}


function escapeHtml(value) {
  return String(value)
    .replaceAll("&", "&amp;")
    .replaceAll("<", "&lt;")
    .replaceAll(">", "&gt;")
    .replaceAll('"', "&quot;")
    .replaceAll("'", "&#39;");
}


initialize().catch((error) => {
  showError(error.message || "Failed to initialize the application.");
  setStatus("Initialization failed.");
});
