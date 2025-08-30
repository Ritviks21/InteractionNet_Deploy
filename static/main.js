// static/main.js

// This function runs as soon as the page loads
$(document).ready(function() {
    
    // Fetch the list of drug names from our new backend endpoint
    $.get('/get_drug_names', function(drugNames) {
        const drug1Select = $('#drug1');
        const drug2Select = $('#drug2');
        // For each drug name, create a new <option> and add it to both dropdowns
        $.each(drugNames, function(index, name) {
            drug1Select.append($('<option>', { value: name, text: name }));
            drug2Select.append($('<option>', a= { value: name, text: name }));
        });
    });

    // The rest of the code is similar, but gets values from the dropdowns
    $('#predict-btn').on('click', function() {
        const drug1Name = $('#drug1').val();
        const drug2Name = $('#drug2').val();

        if (!drug1Name || !drug2Name) {
            showResult('Please select a drug from both dropdowns.', 'error');
            return;
        }

        showResult('Analyzing...', 'loading');
        
        // This function call will render the 3D models
        renderMolecule('viewer1', drug1Name);
        renderMolecule('viewer2', drug2Name);

        // This is where the frontend talks to the backend
        $.ajax({
            url: '/predict',
            type: 'POST',
            contentType: 'application/json',
            data: JSON.stringify({ 'drug1_name': drug1Name, 'drug2_name': drug2Name }),
            success: function(data) {
                if (data.prediction === 1) {
                    showResult('Prediction: High Risk of Interaction', 'error');
                } else {
                    showResult('Prediction: Low Risk of Interaction', 'success');
                }
            },
            error: function(xhr) {
                const errorMessage = xhr.responseJSON ? xhr.responseJSON.error : "Failed to connect to the server.";
                showResult(errorMessage, 'error');
            }
        });
    });

    function renderMolecule(elementId, drugName) {
        let element = $('#' + elementId);
        element.empty();
        let viewer = $3Dmol.createViewer(element, { backgroundColor: 'white' });
        
        // This fetches the 3D data (SDF format) for the interactive model
        $.get(`https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/${drugName}/SDF`, function(sdfData) {
            viewer.addModel(sdfData, 'sdf');
            viewer.setStyle({}, {stick: {}});
            viewer.zoomTo();
            viewer.render();
        }).fail(function() {
            element.html(`<p style="padding: 20px;">Could not load 3D model for "${drugName}"</p>`);
        });
    }

    function showResult(message, type) {
        const resultArea = $('#result-area');
        resultArea.text(message);
        // This logic changes the color of the result box
        resultArea.removeClass('success error loading').addClass(type);
    }
});
