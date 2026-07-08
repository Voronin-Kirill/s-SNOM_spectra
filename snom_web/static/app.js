const state = {
  latestResult: null,
  latestAnalysis: null,
  syncingRange: false,
  currentCalculationId: null,
  stopRequested: false,
  history: [],
  nextSpectrumNumber: 1,
};

const MAX_OSCILLATORS = 10;
const HISTORY_LIMIT = 10;

const materialCatalog = {
  bulk: [],
  layer: [],
  substrate: [],
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

const HISTORY_COLORS = [
  "#3157d8",
  "#d43f6c",
  "#0f8b8d",
  "#7b4cc2",
  "#c77600",
  "#2f855a",
  "#b83280",
  "#475467",
  "#8a5a00",
  "#005f73",
];

const plotIds = ["amplitude-plot", "phase-plot", "permittivity-plot", "residual-plot"];

const form = document.getElementById("calculator-form");
const calculateButton = document.getElementById("calculate-button");
const stopButton = document.getElementById("stop-button");
const downloadResultsButton = document.getElementById("download-results-button");
const formStatus = document.getElementById("form-status");
const errorBox = document.getElementById("error-box");
const builtInMaterialSelect = document.getElementById("built-in-material");
const sampleBuiltIn = document.getElementById("sample-built-in");
const sampleUpload = document.getElementById("sample-upload");
const sampleDrudeLorentz = document.getElementById("sample-drude-lorentz");
const bulkSamplePanel = document.getElementById("bulk-sample-panel");
const layeredSamplePanel = document.getElementById("layered-sample-panel");
const layeredMaterials = document.getElementById("layered-materials");
const layeredStructureFields = document.getElementById("layered-structure-fields");
const modelSchemeImage = document.getElementById("model-scheme-image");
const modelSchemeCaption = document.getElementById("model-scheme-caption");
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
const useDrudeTerm = document.getElementById("dl-use-drude");
const drudeTermFields = document.getElementById("dl-drude-fields");
const oscillatorList = document.getElementById("oscillator-list");
const addOscillatorButton = document.getElementById("add-oscillator-button");
const imageChargeCount = document.getElementById("image-charge-count");
const amplitudeCardTitle = document.getElementById("amplitude-card-title");
const phaseCardTitle = document.getElementById("phase-card-title");

const layeredSlotDefinitions = [
  { key: "layer1", title: "Layer 1", role: "layer", hasThickness: true, defaultThickness: 20 },
  { key: "layer2", title: "Layer 2", role: "layer", hasThickness: true, defaultThickness: 50 },
  { key: "substrate", title: "Substrate", role: "substrate", hasThickness: false },
];


async function initialize() {
  await loadMaterials();
  renderLayeredMaterialPanels();
  addOscillator({ strength: 1, resonanceFrequency: 1700, damping: 20 });
  syncUiState();

  form.addEventListener("change", syncUiState);
  form.addEventListener("submit", handleSubmit);
  stopButton.addEventListener("click", handleStopCalculation);
  downloadResultsButton.addEventListener("click", handleDownloadResults);
  addOscillatorButton.addEventListener("click", () => addOscillator());
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
  materialCatalog.bulk = payload.materials || [];
  materialCatalog.layer = payload.layerMaterials || materialCatalog.bulk;
  materialCatalog.substrate = payload.substrateMaterials || materialCatalog.bulk;

  builtInMaterialSelect.innerHTML = "";
  for (const material of materialCatalog.bulk) {
    const option = document.createElement("option");
    option.value = material.id;
    option.textContent = material.label;
    builtInMaterialSelect.appendChild(option);
  }
}


function renderLayeredMaterialPanels() {
  layeredMaterials.innerHTML = "";
  for (const definition of layeredSlotDefinitions) {
    const panel = document.createElement("section");
    panel.className = "subcard material-editor";
    panel.dataset.slot = definition.key;
    panel.dataset.role = definition.role;
    panel.innerHTML = `
      <div class="subcard-header material-editor-header">
        <h3>${escapeHtml(definition.title)}</h3>
      </div>
      <div class="field-grid two-columns material-thickness-row ${definition.hasThickness ? "" : "hidden"}">
        <label>
          <span>Thickness d${definition.key === "layer1" ? "1" : "2"}, nm</span>
          <input class="material-thickness" type="number" value="${definition.defaultThickness || 20}" step="0.1" min="0">
        </label>
      </div>
      <div class="choice-grid sample-method-grid material-method-grid">
        <label class="choice-card active">
          <input type="radio" name="${definition.key}InputMethod" value="builtIn" checked>
          <span>Built-in material</span>
        </label>
        <label class="choice-card">
          <input type="radio" name="${definition.key}InputMethod" value="upload">
          <span>Upload file</span>
        </label>
        <label class="choice-card">
          <input type="radio" name="${definition.key}InputMethod" value="drudeLorentz">
          <span>Drude-Lorentz model</span>
        </label>
      </div>
      <div class="material-built-in field-grid single-column">
        <label>
          <span>Choose the material</span>
          <select class="material-built-in-select"></select>
        </label>
      </div>
      <div class="material-upload field-grid single-column hidden">
        <label>
          <span>Upload permittivity file</span>
          <input class="material-permittivity-file" type="file" accept=".txt,.csv,.dat">
        </label>
        <p class="help-text">Expected format: <code>frequency_cm-1 Re_epsilon Im_epsilon</code></p>
        <a class="inline-link" href="/api/example/permittivity">Download example file</a>
      </div>
      <div class="material-drude-lorentz drude-lorentz-panel hidden">
        <p class="help-text formula-text">
          ε(ω) = ε<sub>∞</sub> − ω<sub>p</sub><sup>2</sup>/(ω<sup>2</sup> + iγ<sub>D</sub>ω)
          + &sum; f<sub>j</sub>ω<sub>j</sub><sup>2</sup>/(ω<sub>j</sub><sup>2</sup> − ω<sup>2</sup> − iγ<sub>j</sub>ω)
        </p>
        <div class="field-grid single-column">
          <label>
            <span>High-frequency permittivity ε<sub>∞</sub></span>
            <input class="material-dl-epsilon-infinity" type="number" value="1" step="0.1" min="0">
          </label>
        </div>
        <div class="subcard drude-card">
          <label class="inline-toggle">
            <input class="material-dl-use-drude" type="checkbox">
            <span>Use Drude term</span>
          </label>
          <div class="field-grid compact two-columns material-dl-drude-fields is-disabled">
            <label>
              <span>Plasma frequency ω<sub>p</sub>, cm<sup>-1</sup></span>
              <input class="material-dl-plasma-frequency" type="number" value="1000" step="1" min="0">
            </label>
            <label>
              <span>Damping γ<sub>D</sub>, cm<sup>-1</sup></span>
              <input class="material-dl-drude-damping" type="number" value="100" step="1" min="0">
            </label>
          </div>
        </div>
        <div class="subcard oscillator-card">
          <div class="subcard-header">
            <h3>Lorentz oscillators</h3>
            <button type="button" class="secondary compact-button material-add-oscillator">Add oscillator</button>
          </div>
          <div class="oscillator-list material-oscillator-list"></div>
          <p class="help-text">Maximum 10 oscillators. Use Drude only, Lorentz only, or both.</p>
        </div>
      </div>
    `;

    const select = panel.querySelector(".material-built-in-select");
    const catalog = definition.role === "substrate" ? materialCatalog.substrate : materialCatalog.layer;
    for (const material of catalog) {
      const option = document.createElement("option");
      option.value = material.id;
      option.textContent = material.label;
      select.appendChild(option);
    }

    const oscillatorListElement = panel.querySelector(".material-oscillator-list");
    panel.querySelector(".material-add-oscillator").addEventListener("click", () => {
      addOscillatorToList(oscillatorListElement);
    });
    addOscillatorToList(oscillatorListElement, { strength: 1, resonanceFrequency: 950, damping: 20 });
    layeredMaterials.appendChild(panel);
  }
}


function syncUiState() {
  const modelType = getCheckedValue("modelType");
  const sampleInputMethod = getCheckedValue("sampleInputMethod");
  const algorithm = getCheckedValue("algorithm");
  const layeredSelected = modelType === "layered";
  const layeredStructure = getCheckedValue("layeredStructure");

  document.querySelectorAll(".choice-card").forEach((card) => {
    const input = card.querySelector("input[type='radio']");
    card.classList.toggle("active", Boolean(input && input.checked));
  });

  sampleBuiltIn.classList.toggle("hidden", sampleInputMethod !== "builtIn");
  sampleUpload.classList.toggle("hidden", sampleInputMethod !== "upload");
  sampleDrudeLorentz.classList.toggle("hidden", sampleInputMethod !== "drudeLorentz");
  bulkSamplePanel.classList.toggle("hidden", layeredSelected);
  layeredSamplePanel.classList.toggle("hidden", !layeredSelected);
  layeredStructureFields.classList.toggle("hidden", !layeredSelected);

  if (layeredSelected) {
    modelSchemeImage.src = "/assets/tip-scheme-layered.jpg";
    modelSchemeImage.alt = "Layered tip geometry scheme";
    modelSchemeCaption.textContent = "Layered sample geometry scheme.";
  } else {
    modelSchemeImage.src = "/assets/tip-scheme-bulk.jpg";
    modelSchemeImage.alt = "Bulk tip geometry scheme";
    modelSchemeCaption.textContent = "Bulk sample geometry scheme.";
  }

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
  imageChargeCount.disabled = !layeredSelected;

  const drudeEnabled = useDrudeTerm.checked;
  drudeTermFields.classList.toggle("is-disabled", !drudeEnabled);
  drudeTermFields.querySelectorAll("input").forEach((input) => {
    input.disabled = !drudeEnabled;
  });
  addOscillatorButton.disabled = oscillatorList.children.length >= MAX_OSCILLATORS;

  document.querySelectorAll(".material-editor").forEach((panel) => {
    const isLayer2 = panel.dataset.slot === "layer2";
    panel.classList.toggle("hidden", !layeredSelected || (isLayer2 && layeredStructure !== "double"));
    syncMaterialEditor(panel);
  });
}


function syncMaterialEditor(panel) {
  const inputMethod = panel.querySelector("input[type='radio']:checked")?.value || "builtIn";
  panel.querySelector(".material-built-in").classList.toggle("hidden", inputMethod !== "builtIn");
  panel.querySelector(".material-upload").classList.toggle("hidden", inputMethod !== "upload");
  panel.querySelector(".material-drude-lorentz").classList.toggle("hidden", inputMethod !== "drudeLorentz");

  const useDrude = panel.querySelector(".material-dl-use-drude");
  const drudeFields = panel.querySelector(".material-dl-drude-fields");
  drudeFields.classList.toggle("is-disabled", !useDrude.checked);
  drudeFields.querySelectorAll("input").forEach((input) => {
    input.disabled = !useDrude.checked;
  });

  const oscillatorListElement = panel.querySelector(".material-oscillator-list");
  panel.querySelector(".material-add-oscillator").disabled = oscillatorListElement.children.length >= MAX_OSCILLATORS;
}


async function handleSubmit(event) {
  event.preventDefault();
  clearError();
  state.currentCalculationId = createCalculationId();
  state.stopRequested = false;
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
    const result = await response.json().catch(() => ({}));

    if (!response.ok) {
      throw new Error(result.error || "Calculation failed.");
    }

    addResultToHistory(result);
    renderResult(result);
    downloadResultsButton.hidden = false;
    setStatus("Calculation finished.");
  } catch (error) {
    if (state.stopRequested && error.message === "Calculation stopped.") {
      clearError();
      setStatus("Calculation stopped.");
      return;
    }
    showError(error.message || "Calculation failed.");
    setStatus("Calculation failed.");
  } finally {
    setBusy(false, "");
    state.currentCalculationId = null;
    state.stopRequested = false;
  }
}


async function handleStopCalculation() {
  if (!state.currentCalculationId) {
    return;
  }

  state.stopRequested = true;
  stopButton.disabled = true;
  setStatus("Stopping calculation...");

  try {
    await fetch("/api/cancel", {
      method: "POST",
      headers: {
        "Content-Type": "application/json",
      },
      body: JSON.stringify({ calculationId: state.currentCalculationId }),
    });
  } catch (error) {
    showError(error.message || "Failed to stop calculation.");
  }
}


async function buildPayload() {
  const modelType = getCheckedValue("modelType");
  const sampleInputMethod = getCheckedValue("sampleInputMethod");
  const experimentalOn = experimentalEnabled.checked;

  return {
    calculationId: state.currentCalculationId,
    modelType,
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
    sample: modelType === "bulk" ? await buildBulkSamplePayload(sampleInputMethod) : null,
    layered: modelType === "layered" ? await buildLayeredPayload() : null,
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


async function buildBulkSamplePayload(sampleInputMethod) {
  return {
    inputMethod: sampleInputMethod,
    builtInMaterial: sampleInputMethod === "builtIn" ? builtInMaterialSelect.value : null,
    uploadedPermittivity:
      sampleInputMethod === "upload"
        ? await readFileContent(document.getElementById("permittivity-file"))
        : null,
    drudeLorentz: sampleInputMethod === "drudeLorentz" ? collectDrudeLorentzParameters() : null,
  };
}


async function buildLayeredPayload() {
  const structure = getCheckedValue("layeredStructure");
  const layer1Panel = getMaterialPanel("layer1");
  const layer2Panel = getMaterialPanel("layer2");
  const substratePanel = getMaterialPanel("substrate");
  return {
    structure,
    imageCharges: getInteger("image-charge-count"),
    layer1: {
      thickness: getPanelNumber(layer1Panel, ".material-thickness"),
      material: await collectMaterialInput(layer1Panel),
    },
    layer2:
      structure === "double"
        ? {
            thickness: getPanelNumber(layer2Panel, ".material-thickness"),
            material: await collectMaterialInput(layer2Panel),
          }
        : null,
    substrate: await collectMaterialInput(substratePanel),
  };
}


async function collectMaterialInput(panel) {
  const inputMethod = panel.querySelector("input[type='radio']:checked").value;
  return {
    inputMethod,
    builtInMaterial:
      inputMethod === "builtIn" ? panel.querySelector(".material-built-in-select").value : null,
    uploadedPermittivity:
      inputMethod === "upload"
        ? await readFileContent(panel.querySelector(".material-permittivity-file"))
        : null,
    drudeLorentz:
      inputMethod === "drudeLorentz" ? collectMaterialDrudeLorentzParameters(panel) : null,
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
  if (!isPositiveInteger(payload.harmonicNumber)) {
    errors.push("Harmonic number must be a positive integer.");
  }
  if (!(payload.frequencyRange.start < payload.frequencyRange.end)) {
    errors.push("Start frequency must be smaller than end frequency.");
  }
  if (!isPositiveInteger(payload.frequencyRange.points) || payload.frequencyRange.points > 1000) {
    errors.push("Number of frequency points must be between 1 and 1000.");
  }
  if (!(payload.farField.angleDegrees >= 0 && payload.farField.angleDegrees <= 90)) {
    errors.push("Illumination angle must be between 0 and 90 degrees.");
  }

  if (payload.algorithm === "accurate") {
    if (!isPositiveInteger(payload.discretization.N) || payload.discretization.N > 500) {
      errors.push("N must be between 1 and 500.");
    }
    if (!isPositiveInteger(payload.discretization.M) || payload.discretization.M > 100) {
      errors.push("M must be between 1 and 100.");
    }
  }

  if (payload.modelType === "bulk") {
    validateMaterialInput(payload.sample, "Sample", errors);
  }

  if (payload.modelType === "layered") {
    validateLayeredPayload(payload.layered, errors);
  }

  if (payload.experimentalSpectrum.enabled && !payload.experimentalSpectrum.content) {
    errors.push("Please upload an experimental spectrum file.");
  }

  return errors;
}


function validateLayeredPayload(layered, errors) {
  if (!layered) {
    errors.push("Layered sample parameters are required.");
    return;
  }
  if (!(layered.imageCharges > 0 && Number.isInteger(layered.imageCharges) && layered.imageCharges <= 50)) {
    errors.push("K must be between 1 and 50.");
  }
  if (!(layered.layer1.thickness > 0)) {
    errors.push("Layer 1 thickness must be positive.");
  }
  validateMaterialInput(layered.layer1.material, "Layer 1", errors);

  if (layered.structure === "double") {
    if (!(layered.layer2 && layered.layer2.thickness > 0)) {
      errors.push("Layer 2 thickness must be positive.");
    } else {
      validateMaterialInput(layered.layer2.material, "Layer 2", errors);
    }
  }

  validateMaterialInput(layered.substrate, "Substrate", errors);
}


function validateMaterialInput(material, label, errors) {
  if (!material) {
    errors.push(`${label} material is required.`);
    return;
  }
  if (material.inputMethod === "upload" && !material.uploadedPermittivity) {
    errors.push(`Please upload a permittivity file for ${label}.`);
  }
  if (material.inputMethod === "drudeLorentz") {
    validateDrudeLorentzParameters(material.drudeLorentz, errors, label);
  }
}


function validateDrudeLorentzParameters(parameters, errors, label = "Drude-Lorentz") {
  if (!parameters) {
    errors.push(`${label} Drude-Lorentz parameters are required.`);
    return;
  }
  if (!(parameters.epsilonInfinity > 0)) {
    errors.push(`${label}: high-frequency permittivity must be greater than 0.`);
  }

  if (parameters.useDrude) {
    if (!(parameters.plasmaFrequency > 0)) {
      errors.push(`${label}: Drude plasma frequency must be greater than 0.`);
    }
    if (!(parameters.drudeDamping > 0)) {
      errors.push(`${label}: Drude damping must be greater than 0.`);
    }
  }

  if (parameters.oscillators.length > MAX_OSCILLATORS) {
    errors.push(`At most ${MAX_OSCILLATORS} Lorentz oscillators are allowed.`);
  }
  if (!parameters.useDrude && parameters.oscillators.length === 0) {
    errors.push(`${label}: add at least one Lorentz oscillator or enable the Drude term.`);
  }

  parameters.oscillators.forEach((oscillator, index) => {
    const label = `Oscillator ${index + 1}`;
    if (!(oscillator.strength > 0)) {
      errors.push(`${label} strength must be greater than 0.`);
    }
    if (!(oscillator.resonanceFrequency > 0)) {
      errors.push(`${label} resonance frequency must be greater than 0.`);
    }
    if (!(oscillator.damping > 0)) {
      errors.push(`${label} damping must be greater than 0.`);
    }
  });
}


function collectDrudeLorentzParameters() {
  return {
    epsilonInfinity: getNumber("dl-epsilon-infinity"),
    useDrude: useDrudeTerm.checked,
    plasmaFrequency: getNumber("dl-plasma-frequency"),
    drudeDamping: getNumber("dl-drude-damping"),
    oscillators: collectOscillators(oscillatorList),
  };
}


function collectMaterialDrudeLorentzParameters(panel) {
  return {
    epsilonInfinity: getPanelNumber(panel, ".material-dl-epsilon-infinity"),
    useDrude: panel.querySelector(".material-dl-use-drude").checked,
    plasmaFrequency: getPanelNumber(panel, ".material-dl-plasma-frequency"),
    drudeDamping: getPanelNumber(panel, ".material-dl-drude-damping"),
    oscillators: collectOscillators(panel.querySelector(".material-oscillator-list")),
  };
}


function collectOscillators(list) {
  return Array.from(list.querySelectorAll(".oscillator-row")).map((row) => ({
    strength: Number(row.querySelector('[data-field="strength"]').value),
    resonanceFrequency: Number(row.querySelector('[data-field="resonanceFrequency"]').value),
    damping: Number(row.querySelector('[data-field="damping"]').value),
  }));
}


function addOscillator(values = {}) {
  addOscillatorToList(oscillatorList, values);
}


function addOscillatorToList(list, values = {}) {
  if (list.children.length >= MAX_OSCILLATORS) {
    return;
  }

  const row = document.createElement("div");
  row.className = "oscillator-row";
  row.innerHTML = `
    <div class="oscillator-row-header">
      <strong>Oscillator</strong>
      <button type="button" class="secondary compact-button remove-oscillator-button">Remove</button>
    </div>
    <div class="field-grid three-columns">
      <label>
        <span>Strength f<sub>j</sub></span>
        <input data-field="strength" type="number" value="${values.strength ?? 1}" step="0.1" min="0">
      </label>
      <label>
        <span>Resonance ω<sub>j</sub>, cm<sup>-1</sup></span>
        <input data-field="resonanceFrequency" type="number" value="${values.resonanceFrequency ?? 1700}" step="1" min="0">
      </label>
      <label>
        <span>Damping γ<sub>j</sub>, cm<sup>-1</sup></span>
        <input data-field="damping" type="number" value="${values.damping ?? 20}" step="1" min="0">
      </label>
    </div>
  `;

  row.querySelector(".remove-oscillator-button").addEventListener("click", () => {
    row.remove();
    renumberOscillators(list);
    syncUiState();
  });
  list.appendChild(row);
  renumberOscillators(list);
  syncUiState();
}


function renumberOscillators(list = oscillatorList) {
  list.querySelectorAll(".oscillator-row").forEach((row, index) => {
    row.querySelector("strong").textContent = `Oscillator ${index + 1}`;
  });
}


function addResultToHistory(result) {
  result.sessionLabel = `Spectrum ${state.nextSpectrumNumber}`;
  state.nextSpectrumNumber += 1;
  state.history.push(result);
  if (state.history.length > HISTORY_LIMIT) {
    state.history.shift();
  }
  state.latestResult = result;
}


function renderResult(result) {
  const analysis = analyzeResult(result);
  state.latestAnalysis = analysis;
  updateResultLabels(result.metadata.harmonicNumber);
  updateSummary(result.metadata, analysis);
  renderPlots(result, analysis);
  updateDataInspector(analysis.maxAmplitudeFrequency);
}


function updateSummary(metadata, analysis) {
  resultSummary.classList.remove("muted-summary");
  const sigmaHtml = `σ<sub>${escapeHtml(metadata.harmonicNumber)}</sub>`;
  resultSummary.innerHTML = `
    <div>
      <span class="summary-label">|${sigmaHtml}| peak</span>
      <strong>${formatNumber(analysis.maxAmplitude, 4)}</strong>
    </div>
    <div>
      <span class="summary-label">|${sigmaHtml}| peak position</span>
      <strong>${formatNumber(analysis.maxAmplitudeFrequency, 2)} cm⁻¹</strong>
    </div>
    <div>
      <span class="summary-label">arg(${sigmaHtml}) peak</span>
      <strong>${formatNumber(analysis.maxPhase, 4)} rad</strong>
    </div>
    <div>
      <span class="summary-label">arg(${sigmaHtml}) peak position</span>
      <strong>${formatNumber(analysis.maxPhaseFrequency, 2)} cm⁻¹</strong>
    </div>
  `;
}


function updateResultLabels(harmonicNumber) {
  const escapedHarmonic = escapeHtml(harmonicNumber);
  const sigmaText = `σ${toSubscript(harmonicNumber)}`;
  amplitudeCardTitle.innerHTML = `Near-field Amplitude, |σ<sub>${escapedHarmonic}</sub>|`;
  phaseCardTitle.innerHTML = `Near-field Phase, arg(σ<sub>${escapedHarmonic}</sub>)`;
  document.getElementById("amplitude-plot").dataset.chartTitle = `Near-field Amplitude, |${sigmaText}|`;
  document.getElementById("phase-plot").dataset.chartTitle = `Near-field Phase, arg(${sigmaText})`;
}


function renderPlots(result, analysis) {
  const x = result.frequency;
  const frequencyRange = [x[0], x[x.length - 1]];
  const sigmaText = `σ${toSubscript(result.metadata.harmonicNumber)}`;
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

  const amplitudeTraces = [];
  const phaseTraces = [];
  const permittivityTraces = [];

  state.history.forEach((historyResult, historyIndex) => {
    const isLatest = historyResult === result;
    const traceSigmaText = `σ${toSubscript(historyResult.metadata.harmonicNumber)}`;
    const traceName = historyResult.sessionLabel || `Spectrum ${historyIndex + 1}`;
    const traceColor = historyColor(historyIndex);
    const visible = isLatest ? true : "legendonly";

    amplitudeTraces.push({
      x: historyResult.frequency,
      y: historyResult.amplitude,
      name: `${traceName} |${traceSigmaText}|`,
      mode: "lines",
      visible,
      line: { color: traceColor, width: isLatest ? 2.8 : 2.0 },
      hovertemplate: `Frequency: %{x:.3f} cm⁻¹<br>|${traceSigmaText}|: %{y:.6g}<extra>${traceName}</extra>`,
    });

    phaseTraces.push({
      x: historyResult.frequency,
      y: historyResult.phase,
      name: `${traceName} arg(${traceSigmaText})`,
      mode: "lines",
      visible,
      line: { color: traceColor, width: isLatest ? 2.3 : 1.8 },
      hovertemplate: `Frequency: %{x:.3f} cm⁻¹<br>arg(${traceSigmaText}): %{y:.6g} rad<extra>${traceName}</extra>`,
    });

    if (isLatest) {
      amplitudeTraces.push({
        x: [analysis.maxAmplitudeFrequency],
        y: [analysis.maxAmplitude],
        name: `Peak |${traceSigmaText}|`,
        mode: "markers+text",
        marker: { color: COLORS.marker, size: 7, symbol: "diamond" },
        text: [`peak ${formatNumber(analysis.maxAmplitude, 3)}`],
        textposition: "top center",
        textfont: { size: 10, color: COLORS.marker },
        hovertemplate: `Peak<br>Frequency: %{x:.3f} cm⁻¹<br>|${traceSigmaText}|: %{y:.6g}<extra></extra>`,
      });
    }

    const seriesList = getPermittivitySeries(historyResult);
    seriesList.forEach((series, seriesIndex) => {
      const seriesColor = historyColor(historyIndex + seriesIndex + 1);
      permittivityTraces.push({
        x: historyResult.frequency,
        y: series.epsilonReal,
        name: `${traceName} Re(ε ${series.label})`,
        mode: "lines",
        visible,
        line: { color: seriesColor, width: isLatest ? 2.0 : 1.5 },
        hovertemplate: `Frequency: %{x:.3f} cm⁻¹<br>Re(ε): %{y:.6g}<extra>${traceName} ${series.label}</extra>`,
      });
      permittivityTraces.push({
        x: historyResult.frequency,
        y: series.epsilonImag,
        name: `${traceName} Im(ε ${series.label})`,
        mode: "lines",
        visible,
        line: { color: seriesColor, width: isLatest ? 2.0 : 1.5, dash: "dash" },
        hovertemplate: `Frequency: %{x:.3f} cm⁻¹<br>Im(ε): %{y:.6g}<extra>${traceName} ${series.label}</extra>`,
      });
    });
  });

  if (result.experimental) {
    amplitudeTraces.push({
      x: result.experimental.frequency,
      y: result.experimental.amplitude,
      name: `Experiment |${sigmaText}|`,
      mode: "lines+markers",
      marker: { size: 5, color: COLORS.experiment },
      line: { color: COLORS.experiment, width: 1.6, dash: "dash" },
      hovertemplate: `Frequency: %{x:.3f} cm⁻¹<br>|${sigmaText}|: %{y:.6g}<extra>experiment</extra>`,
    });
    phaseTraces.push({
      x: result.experimental.frequency,
      y: result.experimental.phase,
      name: `Experiment arg(${sigmaText})`,
      mode: "lines+markers",
      marker: { size: 5, color: COLORS.experiment },
      line: { color: COLORS.experiment, width: 1.6, dash: "dash" },
      hovertemplate: `Frequency: %{x:.3f} cm⁻¹<br>arg(${sigmaText}): %{y:.6g} rad<extra>experiment</extra>`,
    });
  }

  Plotly.react(
    "amplitude-plot",
    amplitudeTraces,
    baseLayout("Frequency, ω (cm⁻¹)", `Near-field amplitude, |${sigmaText}|`, frequencyRange),
    plotConfig,
  );

  Plotly.react(
    "phase-plot",
    phaseTraces,
    baseLayout("Frequency, ω (cm⁻¹)", `Near-field phase, arg(${sigmaText}) (rad)`, frequencyRange),
    plotConfig,
  );

  Plotly.react(
    "permittivity-plot",
    permittivityTraces,
    baseLayout("Frequency, ω (cm⁻¹)", "Permittivity, ε", frequencyRange),
    plotConfig,
  );

  if (analysis.comparison) {
    const residualLayout = baseLayout("Frequency, ω (cm⁻¹)", "Amplitude residual", frequencyRange);
    residualLayout.yaxis2 = {
      title: "Phase residual (rad)",
      overlaying: "y",
      side: "right",
      titlefont: { size: 10 },
      tickfont: { size: 9 },
      gridcolor: "rgba(0,0,0,0)",
      zerolinecolor: "#e2e8f4",
    };
    Plotly.react(
      "residual-plot",
      [
        {
          x: analysis.comparison.frequency,
          y: analysis.comparison.amplitudeResidual,
          name: "Amplitude residual",
          mode: "lines+markers",
          marker: { color: COLORS.residual, size: 4 },
          line: { color: COLORS.residual, width: 1.7 },
          hovertemplate: "Frequency: %{x:.3f} cm⁻¹<br>Amplitude residual: %{y:.6g}<extra></extra>",
        },
        {
          x: analysis.comparison.frequency,
          y: analysis.comparison.phaseResidual,
          name: "Phase residual",
          mode: "lines+markers",
          yaxis: "y2",
          marker: { color: COLORS.imaginary, size: 4 },
          line: { color: COLORS.imaginary, width: 1.7, dash: "dash" },
          hovertemplate: "Frequency: %{x:.3f} cm⁻¹<br>Phase residual: %{y:.6g} rad<extra></extra>",
        },
      ],
      residualLayout,
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


function baseLayout(xLabel, yLabel, frequencyRange) {
  return {
    paper_bgcolor: "rgba(0,0,0,0)",
    plot_bgcolor: "#ffffff",
    autosize: true,
    margin: { l: 42, r: 10, t: 12, b: 34 },
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
      y: 1.08,
      font: { size: 9 },
    },
  };
}


function historyColor(index) {
  return HISTORY_COLORS[index % HISTORY_COLORS.length];
}


function getPermittivitySeries(result) {
  if (Array.isArray(result.permittivitySeries) && result.permittivitySeries.length > 0) {
    return result.permittivitySeries;
  }
  return [
    {
      key: "sample",
      label: "sample",
      epsilonReal: result.epsilonReal,
      epsilonImag: result.epsilonImag,
    },
  ];
}


function analyzeResult(result) {
  const peakIndex = indexOfMax(result.amplitude);
  const phasePeakIndex = indexOfMax(result.phase);
  const comparison = result.experimental ? computeExperimentalComparison(result) : null;
  return {
    frequencyMin: result.frequency[0],
    frequencyMax: result.frequency[result.frequency.length - 1],
    maxAmplitude: result.amplitude[peakIndex],
    maxAmplitudeFrequency: result.frequency[peakIndex],
    maxPhase: result.phase[phasePeakIndex],
    maxPhaseFrequency: result.frequency[phasePeakIndex],
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
  const simulatedPhase = result.phase;
  const experimentalFrequency = result.experimental.frequency;
  const experimentalAmplitude = result.experimental.amplitude;
  const experimentalPhase = result.experimental.phase;
  const overlapMin = Math.max(simulatedFrequency[0], experimentalFrequency[0]);
  const overlapMax = Math.min(
    simulatedFrequency[simulatedFrequency.length - 1],
    experimentalFrequency[experimentalFrequency.length - 1],
  );

  const frequency = [];
  const amplitudeResidual = [];
  const phaseResidual = [];
  const simulation = [];
  const experiment = [];

  for (let index = 0; index < simulatedFrequency.length; index += 1) {
    const x = simulatedFrequency[index];
    if (x < overlapMin || x > overlapMax) {
      continue;
    }
    const experimentalAmplitudeValue = interpolateLinear(experimentalFrequency, experimentalAmplitude, x);
    const experimentalPhaseValue = interpolateLinear(experimentalFrequency, experimentalPhase, x);
    if (!Number.isFinite(experimentalAmplitudeValue) || !Number.isFinite(experimentalPhaseValue)) {
      continue;
    }
    const simulatedAmplitudeValue = simulatedAmplitude[index];
    const simulatedPhaseValue = simulatedPhase[index];
    frequency.push(x);
    simulation.push(simulatedAmplitudeValue);
    experiment.push(experimentalAmplitudeValue);
    amplitudeResidual.push(simulatedAmplitudeValue - experimentalAmplitudeValue);
    phaseResidual.push(simulatedPhaseValue - experimentalPhaseValue);
  }

  if (amplitudeResidual.length < 2) {
    return null;
  }

  const absResidual = amplitudeResidual.map(Math.abs);
  const rmse = Math.sqrt(
    amplitudeResidual.reduce((sum, value) => sum + value * value, 0) / amplitudeResidual.length,
  );
  const mae = absResidual.reduce((sum, value) => sum + value, 0) / absResidual.length;
  const correlation = pearsonCorrelation(simulation, experiment);
  const simulatedPeakFrequency = simulatedFrequency[indexOfMax(simulatedAmplitude)];
  const experimentalPeakFrequency = experimentalFrequency[indexOfMax(experimentalAmplitude)];

  return {
    frequency,
    residual: amplitudeResidual,
    amplitudeResidual,
    phaseResidual,
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
  const sigmaText = `σ${toSubscript(result.metadata.harmonicNumber)}`;
  const epsilonRows = getPermittivitySeries(result)
    .map(
      (series) => `
        <tr><th>Re(ε ${escapeHtml(series.label)})</th><td>${formatNumber(series.epsilonReal[index], 6)}</td></tr>
        <tr><th>Im(ε ${escapeHtml(series.label)})</th><td>${formatNumber(series.epsilonImag[index], 6)}</td></tr>
      `,
    )
    .join("");
  dataInspectorBody.innerHTML = `
    <tr><th>Frequency</th><td>${formatNumber(result.frequency[index], 3)} cm⁻¹</td></tr>
    ${epsilonRows}
    <tr><th>|${sigmaText}|</th><td>${formatNumber(result.amplitude[index], 6)}</td></tr>
    <tr><th>arg(${sigmaText})</th><td>${formatNumber(result.phase[index], 6)} rad</td></tr>
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
  layout.margin = { l: 78, r: 28, t: 28, b: 68 };
  delete layout.title;
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

  const result = state.latestResult;
  const permittivitySeries = getPermittivitySeries(result);
  const rows = [
    [
      "frequency_cm-1",
      ...permittivitySeries.flatMap((series) => [
        `Re_epsilon_${sanitizeCsvHeader(series.key)}`,
        `Im_epsilon_${sanitizeCsvHeader(series.key)}`,
      ]),
      "Re_sigma_n",
      "Im_sigma_n",
      "abs_sigma_n",
      "arg_sigma_n",
    ],
  ];

  for (let index = 0; index < result.frequency.length; index += 1) {
    rows.push([
      result.frequency[index],
      ...permittivitySeries.flatMap((series) => [
        series.epsilonReal[index],
        series.epsilonImag[index],
      ]),
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
  stopButton.hidden = !isBusy;
  stopButton.disabled = !isBusy || state.stopRequested;
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
  return Number(document.getElementById(id).value);
}


function getMaterialPanel(slot) {
  return layeredMaterials.querySelector(`[data-slot="${slot}"]`);
}


function getPanelNumber(panel, selector) {
  return Number(panel.querySelector(selector).value);
}


function getCheckedValue(name) {
  return document.querySelector(`input[name="${name}"]:checked`).value;
}


function createCalculationId() {
  if (window.crypto && typeof window.crypto.randomUUID === "function") {
    return window.crypto.randomUUID();
  }
  return `${Date.now()}-${Math.random().toString(16).slice(2)}`;
}


function isPositiveInteger(value) {
  return Number.isInteger(value) && value > 0;
}


function toSubscript(value) {
  const digits = {
    "0": "₀",
    "1": "₁",
    "2": "₂",
    "3": "₃",
    "4": "₄",
    "5": "₅",
    "6": "₆",
    "7": "₇",
    "8": "₈",
    "9": "₉",
  };
  return String(value).replace(/[0-9]/g, (digit) => digits[digit]);
}


function sanitizeCsvHeader(value) {
  return String(value).replace(/[^A-Za-z0-9_]+/g, "_").replace(/^_+|_+$/g, "") || "value";
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
