// static/main.js (Final Version with New Output Format)

$(document).ready(function() {
    $('#predict-btn').on('click', function() {
        const drug1Name = $('#drug1').val().trim();
        const drug2Name = $('#drug2').val().trim();

        if (!drug1Name || !drug2Name) {
            showResult('Please enter both drug names.', 'error');
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
                // --- THIS IS THE KEY CHANGE FOR THE NEW FORMAT ---
                const probabilityPercent = (data.probability * 100).toFixed(2); // Format to 2 decimal places
                if (data.prediction === 1) {
                    const message = `Predicted: INTERACTION LIKELY (${probabilityPercent}%)`;
                    showResult(message, 'error');
                } else {
                    const message = `Predicted: NO INTERACTION (${probabilityPercent}%)`;
                    showResult(message, 'success');
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
        resultArea.removeClass('success error loading').addClass(type);
    }
});
