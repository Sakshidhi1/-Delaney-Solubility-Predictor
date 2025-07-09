import streamlit as st
from rdkit import Chem
from rdkit.Chem import Descriptors
import joblib

model = joblib.load("logS_model.pkl")

def featurize(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return [
        Descriptors.MolWt(mol),
        Descriptors.MolLogP(mol),
        Descriptors.NumRotatableBonds(mol),
        Descriptors.TPSA(mol),
        sum(1 for atom in mol.GetAromaticAtoms()) / mol.GetNumAtoms()
    ]

st.title("🧪 Delaney Solubility Predictor")
st.markdown("Enter a SMILES string to predict logS (mol/L)")

smiles = st.text_input("SMILES", value="CCO")

if st.button("Predict"):
    features = featurize(smiles)
    if features:
        prediction = model.predict([features])[0]
        st.success(f"Predicted logS: {round(prediction, 3)}")
    else:
        st.error("Invalid SMILES string.")