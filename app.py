# app.py

from flask import Flask, request, jsonify, render_template
import pandas as pd
import joblib

# Initialize the Flask application
app = Flask(__name__)

# --- Load our simplified model and data once when the server starts ---
print("Loading the AI model and drug data...")
try:
    model = joblib.load('simple_model.pkl')
    drug_data = pd.read_csv('drug_data_for_app.csv')
    # Create a dictionary for fast lookups
    smiles_dict = pd.Series(drug_data.SMILES.values, index=drug_data.Name).to_dict()
    print("Model and data loaded successfully.")
except FileNotFoundError:
    print("ERROR: Could not find 'simple_model.pkl' or 'drug_data_for_app.csv'.")
    exit()

# --- Application Routes ---

# This route serves your index.html file
@app.route('/')
def home():
    return render_template('index.html')

# This is the API endpoint that will make predictions
@app.route('/predict', methods=['POST'])
def predict():
    data = request.get_json()
    drug1_name = data.get('drug1_name')
    drug2_name = data.get('drug2_name')

    # Find the SMILES strings for the given drug names
    smiles1 = smiles_dict.get(drug1_name)
    smiles2 = smiles_dict.get(drug2_name)
        
    if smiles1 is None or smiles2 is None:
        return jsonify({'error': 'One or both drugs not found in our dataset.'}), 400

    # The features for our simple model are the lengths of the SMILES strings
    len1 = len(smiles1)
    len2 = len(smiles2)
    
    # Make a prediction
    # The model expects a 2D array, so we use [[len1, len2]]
    prediction_score = model.predict([[len1, len2]])[0]
    
    # Send the result back to the frontend
    # We'll also send the SMILES strings so the frontend can render them
    return jsonify({
        'prediction_score': round(prediction_score, 2),
        'smiles1': smiles1,
        'smiles2': smiles2
    })

if __name__ == '__main__':
    # This is for local testing
    app.run(debug=True)
