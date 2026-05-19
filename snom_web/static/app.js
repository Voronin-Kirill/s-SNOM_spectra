const state = {
  latestResult: null,
  latestAnalysis: null,
  syncingRange: false,
};

const COLORS = {
  simulation: "#3157d8",
  experiment: "#ff8f3f",
  residual: "#7b8497",
  real: "#3157d8",
  imaginary: "#d43f6c",
  magnitude: "#6f63d9",
  marker: "#101828",
};

const plotIds = ["amplitude-plot", "phase-plot", "permittivity-plot", "residual-plot"];

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
const residualCard = document.getElementById("residual-card");
const resultsPanel = document.querySelector(".results-panel");
const dataInspectorBody = document.getElementById("data-inspector-body");


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
  const analysis = analyzeResult(result);
  state.latestAnalysis = analysis;
  updateSummary(result.metadata, analysis);
  renderPlots(result, analysis);
  updateDataInspector(analysis.maxAmplitudeFrequency);
}


function updateSummary(metadata, analysis) {
  resultSummary.classList.remove("muted-summary");
  const comparison = analysis.comparison;
  resultSummary.innerHTML = `
    <div>
      <span class="summary-label">Material / model</span>
      <strong>${escapeHtml(metadata.materialLabel)} / bulk</strong>
    </div>
    <div>
      <span class="summary-label">Frequency range</span>
      <strong>${formatNumber(analysis.frequencyMin, 0)}–${formatNumber(analysis.frequencyMax, 0)} cm⁻¹</strong>
    </div>
    <div>
      <span class="summary-label">Harmonic</span>
      <strong>n = ${metadata.harmonicNumber}</strong>
    </div>
    <div>
      <span class="summary-label">Max |σₙ|</span>
      <strong>${formatNumber(analysis.maxAmplitude, 4)} @ ${formatNumber(analysis.maxAmplitudeFrequency, 1)} cm⁻¹</strong>
    </div>
    <div>
      <span class="summary-label">Phase range</span>
      <strong>${formatNumber(analysis.phaseMin, 3)} to ${formatNumber(analysis.phaseMax, 3)} rad</strong>
    </div>
    <div>
      <span class="summary-label">Re(ε) range</span>
      <strong>${formatNumber(analysis.epsilonRealMin, 3)} to ${formatNumber(analysis.epsilonRealMax, 3)}</strong>
    </div>
    <div>
      <span class="summary-label">Im(ε) range</span>
      <strong>${formatNumber(analysis.epsilonImagMin, 3)} to ${formatNumber(analysis.epsilonImagMax, 3)}</strong>
    </div>
    <div>
      <span class="summary-label">Experiment fit</span>
      <strong>${comparison ? `RMSE ${formatNumber(comparison.rmse, 4)}, MAE ${formatNumber(comparison.mae, 4)}, r ${formatNumber(comparison.correlation, 3)}, Δpeak ${formatNumber(comparison.peakShift, 2)} cm⁻¹` : "No experiment"}</strong>
    </div>
  `;
}


function renderPlots(result, analysis) {
  const x = result.frequency;
  const frequencyRange = [x[0], x[x.length - 1]];
  const plotConfig = {
    responsive: true,
    displaylogo: false,
    toImageButtonOptions: {
      format: "png",
      filename: "s-snom-spectrum",
      scale: 2,
    },
  };

  residualCard.classList.toggle("hidden", !analysis.comparison);
  resultsPanel.classList.toggle("has-residual", Boolean(analysis.comparison));

  const peakTrace = {
    x: [analysis.maxAmplitudeFrequency],
    y: [analysis.maxAmplitude],
    name: "Peak |σₙ|",
    mode: "markers+text",
    marker: { color: COLORS.marker, size: 7, symbol: "diamond" },
    text: [`peak ${formatNumber(analysis.maxAmplitude, 3)}`],
    textposition: "top center",
    textfont: { size: 10, color: COLORS.marker },
    hovertemplate: "Peak<br>Frequency: %{x:.3f} cm⁻¹<br>|σₙ|: %{y:.6g}<extra></extra>",
  };
  const amplitudeTraces = [
    {
      x,
      y: result.amplitude,
      name: "Simulation |σₙ|",
      mode: "lines",
      line: { color: COLORS.simulation, width: 2.8 },
      hovertemplate: "Frequency: %{x:.3f} cm⁻¹<br>|σₙ|: %{y:.6g}<extra>simulation</extra>",
    },
    peakTrace,
  ];
  const phaseTraces = [
    {
      x,
      y: result.phase,
      name: "Simulation arg(σₙ)",
      mode: "lines",
      line: { color: COLORS.simulation, width: 2.3 },
      hovertemplate: "Frequency: %{x:.3f} cm⁻¹<br>arg(σₙ): %{y:.6g} rad<extra>simulation</extra>",
    },
  ];

  if (result.experimental) {
    amplitudeTraces.push({
      x: result.experimental.frequency,
      y: result.experimental.amplitude,
      name: "Experiment |σₙ|",
      mode: "lines+markers",
      marker: { size: 5, color: COLORS.experiment },
      line: { color: COLORS.experiment, width: 1.6, dash: "dash" },
      hovertemplate: "Frequency: %{x:.3f} cm⁻¹<br>|σₙ|: %{y:.6g}<extra>experiment</extra>",
    });
    phaseTraces.push({
      x: result.experimental.frequency,
      y: result.experimental.phase,
      name: "Experiment arg(σₙ)",
      mode: "lines+markers",
      marker: { size: 5, color: COLORS.experiment },
      line: { color: COLORS.experiment, width: 1.6, dash: "dash" },
      hovertemplate: "Frequency: %{x:.3f} cm⁻¹<br>arg(σₙ): %{y:.6g} rad<extra>experiment</extra>",
    });
  }

  Plotly.react(
    "amplitude-plot",
    amplitudeTraces,
    baseLayout("Near-field amplitude, |σₙ|", "Frequency, ω (cm⁻¹)", "Near-field amplitude, |σₙ|", frequencyRange),
    plotConfig,
  );

  Plotly.react(
    "phase-plot",
    phaseTraces,
    baseLayout("Near-field phase, arg(σₙ)", "Frequency, ω (cm⁻¹)", "Near-field phase, arg(σₙ) (rad)", frequencyRange),
    plotConfig,
  );

  const epsilonMagnitude = result.epsilonReal.map((real, index) => Math.hypot(real, result.epsilonImag[index]));
  Plotly.react(
    "permittivity-plot",
    [
      {
        x,
        y: result.epsilonReal,
        name: "Re(ε)",
        mode: "lines",
        line: { color: COLORS.real, width: 2.2 },
        hovertemplate: "Frequency: %{x:.3f} cm⁻¹<br>Re(ε): %{y:.6g}<extra></extra>",
      },
      {
        x,
        y: result.epsilonImag,
        name: "Im(ε)",
        mode: "lines",
        line: { color: COLORS.imaginary, width: 2.2 },
        hovertemplate: "Frequency: %{x:.3f} cm⁻¹<br>Im(ε): %{y:.6g}<extra></extra>",
      },
      {
        x,
        y: epsilonMagnitude,
        name: "|ε|",
        mode: "lines",
        visible: "legendonly",
        line: { color: COLORS.magnitude, width: 1.7, dash: "dot" },
        hovertemplate: "Frequency: %{x:.3f} cm⁻¹<br>|ε|: %{y:.6g}<extra></extra>",
      },
    ],
    baseLayout("Complex permittivity, ε", "Frequency, ω (cm⁻¹)", "Permittivity, ε", frequencyRange),
    plotConfig,
  );

  if (analysis.comparison) {
    Plotly.react(
      "residual-plot",
      [
        {
          x: analysis.comparison.frequency,
          y: analysis.comparison.residual,
          name: "Simulation − experiment",
          mode: "lines+markers",
          marker: { color: COLORS.residual, size: 4 },
          line: { color: COLORS.residual, width: 1.7 },
          hovertemplate: "Frequency: %{x:.3f} cm⁻¹<br>Residual: %{y:.6g}<extra></extra>",
        },
      ],
      baseLayout("Amplitude residual", "Frequency, ω (cm⁻¹)", "Simulation − experiment", frequencyRange),
      plotConfig,
    );
  } else {
    Plotly.purge(document.getElementById("residual-plot"));
  }

  plotIds.forEach(attachPlotInteractions);
}


function renderEmptyPlots() {
  ["amplitude-plot", "phase-plot", "permittivity-plot"].forEach((plotId) => {
    const plot = document.getElementById(plotId);
    Plotly.purge(plot);
    plot.innerHTML = `
      <div class="plot-empty">
        <div>
          <strong>No calculation yet</strong>
          <span>Set parameters and click Calculate</span>
        </div>
      </div>
    `;
  });
  residualCard.classList.add("hidden");
  resultsPanel.classList.remove("has-residual");
}


function baseLayout(title, xLabel, yLabel, frequencyRange) {
  return {
    title: { text: title, font: { size: 12 } },
    paper_bgcolor: "rgba(0,0,0,0)",
    plot_bgcolor: "#ffffff",
    autosize: true,
    margin: { l: 42, r: 10, t: 30, b: 34 },
    hovermode: "x unified",
    spikedistance: -1,
    xaxis: {
      title: xLabel,
      range: frequencyRange,
      titlefont: { size: 10 },
      tickfont: { size: 9 },
      gridcolor: "#e2e8f4",
      zerolinecolor: "#e2e8f4",
      showspikes: true,
      spikemode: "across",
      spikesnap: "cursor",
      spikecolor: "rgba(16, 24, 40, 0.45)",
      spikethickness: 1,
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


function analyzeResult(result) {
  const peakIndex = indexOfMax(result.amplitude);
  const comparison = result.experimental ? computeExperimentalComparison(result) : null;
  return {
    frequencyMin: result.frequency[0],
    frequencyMax: result.frequency[result.frequency.length - 1],
    maxAmplitude: result.amplitude[peakIndex],
    maxAmplitudeFrequency: result.frequency[peakIndex],
    phaseMin: Math.min(...result.phase),
    phaseMax: Math.max(...result.phase),
    epsilonRealMin: Math.min(...result.epsilonReal),
    epsilonRealMax: Math.max(...result.epsilonReal),
    epsilonImagMin: Math.min(...result.epsilonImag),
    epsilonImagMax: Math.max(...result.epsilonImag),
    comparison,
  };
}


function computeExperimentalComparison(result) {
  const simulatedFrequency = result.frequency;
  const simulatedAmplitude = result.amplitude;
  const experimentalFrequency = result.experimental.frequency;
  const experimentalAmplitude = result.experimental.amplitude;
  const overlapMin = Math.max(simulatedFrequency[0], experimentalFrequency[0]);
  const overlapMax = Math.min(
    simulatedFrequency[simulatedFrequency.length - 1],
    experimentalFrequency[experimentalFrequency.length - 1],
  );

  const frequency = [];
  const residual = [];
  const simulation = [];
  const experiment = [];

  for (let index = 0; index < simulatedFrequency.length; index += 1) {
    const x = simulatedFrequency[index];
    if (x < overlapMin || x > overlapMax) {
      continue;
    }
    const experimentalValue = interpolateLinear(experimentalFrequency, experimentalAmplitude, x);
    if (!Number.isFinite(experimentalValue)) {
      continue;
    }
    const simulatedValue = simulatedAmplitude[index];
    frequency.push(x);
    simulation.push(simulatedValue);
    experiment.push(experimentalValue);
    residual.push(simulatedValue - experimentalValue);
  }

  if (residual.length < 2) {
    return null;
  }

  const absResidual = residual.map(Math.abs);
  const rmse = Math.sqrt(residual.reduce((sum, value) => sum + value * value, 0) / residual.length);
  const mae = absResidual.reduce((sum, value) => sum + value, 0) / absResidual.length;
  const correlation = pearsonCorrelation(simulation, experiment);
  const simulatedPeakFrequency = simulatedFrequency[indexOfMax(simulatedAmplitude)];
  const experimentalPeakFrequency = experimentalFrequency[indexOfMax(experimentalAmplitude)];

  return {
    frequency,
    residual,
    rmse,
    mae,
    correlation,
    peakShift: simulatedPeakFrequency - experimentalPeakFrequency,
  };
}


function interpolateLinear(xValues, yValues, targetX) {
  if (targetX < xValues[0] || targetX > xValues[xValues.length - 1]) {
    return Number.NaN;
  }
  for (let index = 0; index < xValues.length - 1; index += 1) {
    const x0 = xValues[index];
    const x1 = xValues[index + 1];
    if (targetX >= x0 && targetX <= x1) {
      if (x1 === x0) {
        return yValues[index];
      }
      const t = (targetX - x0) / (x1 - x0);
      return yValues[index] + t * (yValues[index + 1] - yValues[index]);
    }
  }
  return Number.NaN;
}


function pearsonCorrelation(seriesA, seriesB) {
  const meanA = mean(seriesA);
  const meanB = mean(seriesB);
  let numerator = 0;
  let denominatorA = 0;
  let denominatorB = 0;
  for (let index = 0; index < seriesA.length; index += 1) {
    const deltaA = seriesA[index] - meanA;
    const deltaB = seriesB[index] - meanB;
    numerator += deltaA * deltaB;
    denominatorA += deltaA * deltaA;
    denominatorB += deltaB * deltaB;
  }
  const denominator = Math.sqrt(denominatorA * denominatorB);
  return denominator === 0 ? Number.NaN : numerator / denominator;
}


function mean(values) {
  return values.reduce((sum, value) => sum + value, 0) / values.length;
}


function indexOfMax(values) {
  let maxIndex = 0;
  for (let index = 1; index < values.length; index += 1) {
    if (values[index] > values[maxIndex]) {
      maxIndex = index;
    }
  }
  return maxIndex;
}


function formatNumber(value, digits) {
  if (!Number.isFinite(value)) {
    return "n/a";
  }
  return Number(value).toLocaleString("en-US", {
    maximumFractionDigits: digits,
    minimumFractionDigits: 0,
  });
}


function attachPlotInteractions(plotId) {
  const plot = document.getElementById(plotId);
  if (!plot || !plot.data || plot.dataset.interactionsAttached === "true") {
    return;
  }
  plot.dataset.interactionsAttached = "true";
  plot.on("plotly_relayout", (eventData) => syncXAxisRange(plotId, eventData));
  plot.on("plotly_hover", (eventData) => {
    const point = eventData.points && eventData.points[0];
    if (point) {
      updateCrosshair(point.x);
    }
  });
  plot.on("plotly_unhover", clearCrosshair);
}


function syncXAxisRange(sourcePlotId, eventData) {
  if (state.syncingRange) {
    return;
  }
  const range = extractXAxisRange(eventData);
  if (!range) {
    return;
  }

  state.syncingRange = true;
  Promise.all(
    getVisiblePlotIds()
      .filter((plotId) => plotId !== sourcePlotId)
      .map((plotId) => Plotly.relayout(document.getElementById(plotId), { "xaxis.range": range })),
  ).finally(() => {
    state.syncingRange = false;
  });
}


function extractXAxisRange(eventData) {
  if (eventData["xaxis.range[0]"] !== undefined && eventData["xaxis.range[1]"] !== undefined) {
    return [eventData["xaxis.range[0]"], eventData["xaxis.range[1]"]];
  }
  if (Array.isArray(eventData["xaxis.range"])) {
    return eventData["xaxis.range"];
  }
  if (eventData["xaxis.autorange"] && state.latestResult) {
    const x = state.latestResult.frequency;
    return [x[0], x[x.length - 1]];
  }
  return null;
}


function updateCrosshair(frequency) {
  updateDataInspector(frequency);
  const shape = {
    type: "line",
    xref: "x",
    yref: "paper",
    x0: frequency,
    x1: frequency,
    y0: 0,
    y1: 1,
    line: { color: "rgba(16, 24, 40, 0.45)", width: 1, dash: "dot" },
  };
  getVisiblePlotIds().forEach((plotId) => {
    Plotly.relayout(document.getElementById(plotId), { shapes: [shape] });
  });
}


function updateDataInspector(frequency) {
  const result = state.latestResult;
  if (!result) {
    return;
  }
  const index = nearestIndex(result.frequency, frequency);
  dataInspectorBody.innerHTML = `
    <tr><th>Frequency</th><td>${formatNumber(result.frequency[index], 3)} cm⁻¹</td></tr>
    <tr><th>Re(ε)</th><td>${formatNumber(result.epsilonReal[index], 6)}</td></tr>
    <tr><th>Im(ε)</th><td>${formatNumber(result.epsilonImag[index], 6)}</td></tr>
    <tr><th>|σₙ|</th><td>${formatNumber(result.amplitude[index], 6)}</td></tr>
    <tr><th>arg(σₙ)</th><td>${formatNumber(result.phase[index], 6)} rad</td></tr>
  `;
}


function nearestIndex(values, target) {
  let bestIndex = 0;
  let bestDistance = Math.abs(values[0] - target);
  for (let index = 1; index < values.length; index += 1) {
    const distance = Math.abs(values[index] - target);
    if (distance < bestDistance) {
      bestIndex = index;
      bestDistance = distance;
    }
  }
  return bestIndex;
}


function clearCrosshair() {
  getVisiblePlotIds().forEach((plotId) => {
    Plotly.relayout(document.getElementById(plotId), { shapes: [] });
  });
}


function getVisiblePlotIds() {
  return plotIds.filter((plotId) => {
    const plot = document.getElementById(plotId);
    return plot && plot.data && !plot.closest(".plot-card")?.classList.contains("hidden");
  });
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
