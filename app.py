# app.py (Final Version with Probability)

from flask import Flask, request, jsonify, render_template
import pandas as pd
import joblib
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

app = Flask(__name__)

# --- Load Model and Data ---
model = joblib.load('ddi_fingerprint_model.pkl')
drug_data = pd.read_csv('ram_friendly_model_data.csv')
print("Model and data loaded successfully.")

# --- Helper Functions ---
def get_morgan_fingerprint(smiles_string, nBits=2048):
    try:
        mol = Chem.MolFromSmiles(smiles_string)
        if mol:
            return np.array(list(AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=nBits)))
        return np.zeros(nBits)
    except:
        return np.zeros(nBits)

def find_smiles(drug_name, data_df):
    drug_name_lower = drug_name.lower()
    row = data_df[data_df['drug_1_name'].str.lower() == drug_name_lower]
    if not row.empty:
        return row.iloc[0]['smiles_1']
    row = data_df[data_df['drug_2_name'].str.lower() == drug_name_lower]
    if not row.empty:
        return row.iloc[0]['smiles_2']
    return None

# --- Application Routes ---
@app.route('/')
def home():
    return render_template('index.html')

@app.route('/predict', methods=['POST'])
def predict():
    data = request.get_json()
    drug1_name = data.get('drug1_name')
    drug2_name = data.get('drug2_name')

    smiles1 = find_smiles(drug1_name, drug_data)
    smiles2 = find_smiles(drug2_name, drug_data)
        
    if smiles1 is None or smiles2 is None:
        return jsonify({'error': 'One or both drugs not found in our dataset.'}), 400

    fp1 = get_morgan_fingerprint(smiles1)
    fp2 = get_morgan_fingerprint(smiles2)
    
    combined_fp = np.concatenate((fp1, fp2)).reshape(1, -1)
    
    # --- THIS IS THE KEY CHANGE: GETTING THE PROBABILITY ---
    prediction = model.predict(combined_fp)[0]
    probability = model.predict_proba(combined_fp)[0]
    
    # The probability of "High Risk" (class 1) is the second value
    high_risk_prob = probability[1]
    
    return jsonify({
        'prediction': int(prediction),
        'probability': float(high_risk_prob) # Sending the probability to the frontend
    })
