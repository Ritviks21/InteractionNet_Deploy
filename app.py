# app.py (Final, Optimized Version)

from flask import Flask, request, jsonify, render_template
import pandas as pd
import joblib
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
import os

app = Flask(__name__)

# --- Global variables to hold our model and data ---
model = None
smiles_dict = {}
drug_list = []

# --- Helper Functions ---
def get_morgan_fingerprint(smiles_string, nBits=2048):
    try:
        mol = Chem.MolFromSmiles(smiles_string)
        if mol:
            return np.array(list(AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=nBits)))
        return np.zeros(nBits)
    except:
        return np.zeros(nBits)

def load_model_and_data():
    """Load the model and prepare the data dictionary. This runs once."""
    global model, smiles_dict, drug_list
    
    print("Loading the AI model and drug data...")
    try:
        model = joblib.load('ddi_fingerprint_model.pkl')
        drug_data = pd.read_csv('ram_friendly_model_data.csv')
        
        # Create a fast-lookup dictionary for SMILES strings
        for _, row in drug_data.iterrows():
            smiles_dict[row['drug_1_name'].lower()] = row['smiles_1']
            smiles_dict[row['drug_2_name'].lower()] = row['smiles_2']
            
        # Create a sorted list of all unique drug names for the dropdowns
        drug_names_1 = pd.Series(drug_data['drug_1_name'].unique())
        drug_names_2 = pd.Series(drug_data['drug_2_name'].unique())
        all_drug_names = pd.concat([drug_names_1, drug_names_2]).unique()
        all_drug_names.sort()
        drug_list = list(all_drug_names)

        print("Model and data loaded successfully.")
    except FileNotFoundError:
        print("ERROR: Could not find model or data files.")
        exit()

# --- Application Routes ---
@app.route('/')
def home():
    return render_template('index.html')

@app.route('/get_drug_names')
def get_drug_names():
    return jsonify(drug_list)

@app.route('/predict', methods=['POST'])
def predict():
    data = request.get_json()
    drug1_name = data.get('drug1_name', '').lower()
    drug2_name = data.get('drug2_name', '').lower()

    # Use the fast dictionary to look up SMILES strings
    smiles1 = smiles_dict.get(drug1_name)
    smiles2 = smiles_dict.get(drug2_name)
        
    if smiles1 is None or smiles2 is None:
        return jsonify({'error': 'One or both drugs not found in our dataset.'}), 400

    fp1 = get_morgan_fingerprint(smiles1)
    fp2 = get_morgan_fingerprint(smiles2)
    
    combined_fp = np.concatenate((fp1, fp2)).reshape(1, -1)
    prediction = model.predict(combined_fp)[0]
    
    return jsonify({'prediction': int(prediction)})

# --- Main execution ---
if __name__ == '__main__':
    # This runs only when you test locally
    load_model_and_data()
    app.run(debug=True)
else:
    # This runs when deployed on Render
    load_model_and_data()
