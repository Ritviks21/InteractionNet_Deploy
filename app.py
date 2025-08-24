# app.py (Final Version for Deployment)

from flask import Flask, request, jsonify, render_template
import pandas as pd
import joblib
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

# Initialize the Flask application
app = Flask(__name__)

# --- Load the correct model and data files from your repository ---
print("Loading the AI model and drug data...")
try:
    # Use the exact filenames from your GitHub repository
    model = joblib.load('ddi_fingerprint_model.pkl')
    drug_data = pd.read_csv('ram_friendly_model_data.csv')
    print("Model and data loaded successfully.")
except FileNotFoundError:
    print("ERROR: Could not find model or data files. Make sure they are in the root of your repository.")
    exit()

# --- Helper Functions for the Fingerprint Model ---

def get_morgan_fingerprint(smiles_string, nBits=2048):
    """Generates a molecular fingerprint from a SMILES string."""
    try:
        mol = Chem.MolFromSmiles(smiles_string)
        if mol is not None:
            return np.array(list(AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=nBits)))
        return np.zeros(nBits)
    except:
        return np.zeros(nBits)

def find_smiles(drug_name, data_df):
    """Finds the SMILES string for a given drug name in the dataset."""
    drug_name_lower = drug_name.lower()
    # Search in the first drug name column
    smiles_row = data_df[data_df['drug_1_name'].str.lower() == drug_name_lower]
    if not smiles_row.empty:
        return smiles_row.iloc[0]['smiles_1']
    # If not found, search in the second drug name column
    smiles_row = data_df[data_df['drug_2_name'].str.lower() == drug_name_lower]
    if not smiles_row.empty:
        return smiles_row.iloc[0]['smiles_2']
    return None

# --- Application Routes ---

@app.route('/')
def home():
    """Serves the main index.html file."""
    return render_template('index.html') , render_template('main.js') 


@app.route('/predict', methods=['POST'])
def predict():
    """Handles the prediction request from the frontend."""
    data = request.get_json()
    drug1_name = data.get('drug1_name')
    drug2_name = data.get('drug2_name')

    smiles1 = find_smiles(drug1_name, drug_data)
    smiles2 = find_smiles(drug2_name, drug_data)
        
    if smiles1 is None or smiles2 is None:
        return jsonify({'error': 'One or both drugs not found in our dataset.'}), 400

    # Generate fingerprints and combine them for the model
    fp1 = get_morgan_fingerprint(smiles1)
    fp2 = get_morgan_fingerprint(smiles2)
    
    combined_fp = np.concatenate((fp1, fp2)).reshape(1, -1)
    
    # Make a prediction (0 for no interaction, 1 for interaction)
    prediction = model.predict(combined_fp)[0]
    
    # Send the result back to the frontend
    return jsonify({'prediction': int(prediction)}) 


if __name__ == '__main__':
    # This part is for local testing and won't be used by Render
    app.run(debug=True)
