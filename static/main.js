// static/main.js (Final Version with Probability)

document.getElementById('predict-btn').addEventListener('click', () => {
    const drug1Name = document.getElementById('drug1').value.trim();
    const drug2Name = document.getElementById('drug2').value.trim();
    
    if (!drug1Name || !drug2Name) {
        showResult('Please enter both drug names.', 'error', null);
        return;
    }

    showResult('Analyzing...', 'loading', null);
    renderMolecule('viewer1', drug1Name); // This function renders the 3D "output models"
    renderMolecule('viewer2', drug2Name); // This function renders the 3D "output models"
    
    fetch('/predict', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ drug1_name: drug1Name, drug2_name: drug2Name })
    })
    .then(response => {
        if (!response.ok) {
            return response.json().then(err => { throw new Error(err.error) });
        }
        return response.json();
    })
    .then(data => {
        if (data.error) {
            showResult(data.error, 'error', null);
        } else {
            const probabilityPercent = (data.probability * 100).toFixed(1);
            if (data.prediction === 1) {
                showResult(`Prediction: High Risk of Interaction`, 'error', probabilityPercent);
            } else {
                // For low risk, we show 100 - probability
                const lowRiskPercent = (100 - (data.probability * 100)).toFixed(1);
                showResult(`Prediction: Low Risk of Interaction`, 'success', lowRiskPercent);
            }
        }
    })
    .catch(error => {
        console.error('Error:', error);
        showResult(error.message, 'error', null);
    });
});

function renderMolecule(elementId, drugName) {
    let element = document.getElementById(elementId);
    element.innerHTML = '';
    let viewer = $3Dmol.createViewer(element, { backgroundColor: 'white' });
    
    fetch(`https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/${drugName}/SDF`)
        .then(res => res.text())
        .then(sdfData => {
            if (sdfData.includes("PUGREST.NotFound")) {
                 element.innerHTML = `<p style="padding: 20px;">3D model not found for "${drugName}"</p>`;
                 return;
            }
            viewer.addModel(sdfData, 'sdf');
            viewer.setStyle({}, {stick: {}});
            viewer.zoomTo();
            viewer.render();
        })
        .catch(err => {
            console.error("Could not load 3D model for", drugName);
            element.innerHTML = `<p style="padding: 20px;">Could not load 3D model.</p>`;
        });
}

function showResult(message, type, probability) {
    const resultArea = document.getElementById('result-area');
    
    // --- THIS IS THE KEY CHANGE: DISPLAYING THE PROBABILITY ---
    let displayText = message;
    if (probability !== null) {
        displayText += ` (Confidence: ${probability}%)`;
    }
    resultArea.textContent = displayText;
    
    resultArea.style.backgroundColor = '';
    resultArea.style.color = '';
    resultArea.style.borderColor = '';

    if (type === 'error') {
        resultArea.style.backgroundColor = '#f8d7da';
        resultArea.style.color = '#721c24';
        resultArea.style.borderColor = '#f5c6cb';
    } else if (type === 'success') {
        resultArea.style.backgroundColor = '#d1e7dd';
        resultArea.style.color = '#0f5132';
        resultArea.style.borderColor = '#badbcc';
    } else { // 'loading'
        resultArea.style.backgroundColor = '#e2e3e5';
        resultArea.style.color = '#41464b';
        resultArea.style.borderColor = '#d6d8db';
    }
}
