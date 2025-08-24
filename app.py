from flask import Flask, request, jsonify, render_template
import pandas as pd
import joblib
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

# Initialize the Flask application
app = Flask(__name__)

# --- Load the model and drug data ---
print("Loading the AI model and drug data...")
try:
    model = joblib.load('ddi_fingerprint_model.pkl')
    drug_data = pd.read_csv('static/ram_friendly_model_data.csv')  # ✅ moved to static folder
    print("Model and data loaded successfully.")
except FileNotFoundError:
    print("ERROR: Could not find model or data files. Make sure they are in the correct location.")
    exit()

# --- Helper Functions ---

def get_morgan_fingerprint(smiles_string, nBits=2048):
    """Generates a molecular fingerprint from a SMILES string."""
    try:
        mol = Chem.MolFromSmiles(smiles_string)
        if mol is not None:
            return np.array(list(AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=nBits)))
        print(f"Warning: Invalid SMILES string '{smiles_string}'")
        return np.zeros(nBits)
    except Exception as e:
        print(f"Error generating fingerprint: {e}")
        return np.zeros(nBits)

def find_smiles(drug_name, data_df):
    """Finds the SMILES string for a given drug name in the dataset."""
    drug_name_lower = drug_name.lower()
    smiles_row = data_df[data_df['drug_1_name'].str.lower() == drug_name_lower]
    if not smiles_row.empty:
        return smiles_row.iloc[0]['smiles_1']
    smiles_row = data_df[data_df['drug_2_name'].str.lower() == drug_name_lower]
    if not smiles_row.empty:
        return smiles_row.iloc[0]['smiles_2']
    return None

# --- Routes ---

@app.route('/')
def home():
    """Serves the main index.html file."""
    return render_template('index.html')

@app.route('/predict', methods=['POST'])
def predict():
    """Handles the prediction request from the frontend."""
    data = request.get_json()
    drug1_name = data.get('drug1_name')
    drug2_name = data.get('drug2_name')

    print(f"Received prediction request for: {drug1_name} and {drug2_name}")

    smiles1 = find_smiles(drug1_name, drug_data)
    smiles2 = find_smiles(drug2_name, drug_data)

    missing = []
    if smiles1 is None:
        missing.append(drug1_name)
    if smiles2 is None:
        missing.append(drug2_name)

    if missing:
        return jsonify({'error': f"Drug(s) not found: {', '.join(missing)}"}), 400

    fp1 = get_morgan_fingerprint(smiles1)
    fp2 = get_morgan_fingerprint(smiles2)
    combined_fp = np.concatenate((fp1, fp2)).reshape(1, -1)

    prediction = model.predict(combined_fp)[0]
    return jsonify({'prediction': int(prediction)})

@app.route('/health')
def health():
    """Simple health check route."""
    return jsonify({'status': 'ok'})

if __name__ == '__main__':
    app.run(debug=True)
