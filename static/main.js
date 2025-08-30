// static/main.js

// This function runs as soon as the page document is ready
$(document).ready(function() {
    
    // Fetch the list of drug names from our backend to populate dropdowns
    $.get('/get_drug_names', function(drugNames) {
        const drug1Select = $('#drug1');
        const drug2Select = $('#drug2');

        // Clear the "Loading..." message
        drug1Select.empty().append($('<option>', { value: "", text: "Select Drug 1..." }));
        drug2Select.empty().append($('<option>', { value: "", text: "Select Drug 2..." }));
        
        // Add each drug name to both dropdowns
        $.each(drugNames, function(index, name) {
            drug1Select.append($('<option>', { value: name, text: name }));
            drug2Select.append($('<option>', { value: name, text: name }));
        });
    }).fail(function() {
        console.error("Failed to load drug names from the server.");
        $('#drug1').empty().append($('<option>', { value: "", text: "Error loading drugs" }));
        $('#drug2').empty().append($('<option>', { value: "", text: "Error loading drugs" }));
    });

    // Set up the click event for the predict button
    $('#predict-btn').on('click', function() {
        const drug1Name = $('#drug1').val();
        const drug2Name = $('#drug2').val();

        if (!drug1Name || !drug2Name) {
            showResult('Please select a drug from both dropdowns.', 'error');
            return;
        }

        showResult('Analyzing...', 'loading');
        renderMolecule('viewer1', drug1Name);
        renderMolecule('viewer2', drug2Name);

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
});

function renderMolecule(elementId, drugName) {
    let element = $('#' + elementId);
    element.empty();
    let viewer = $3Dmol.createViewer(element, { backgroundColor: 'white' });
    
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
    resultArea.attr('class', '').addClass(type); // A better way to set the class
}
