import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
import joblib
from math import sqrt

# Load dataset
df = pd.read_csv("delaney.csv")

# Feature extraction
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

features = []
targets = []

for _, row in df.iterrows():
    feat = featurize(row["SMILES"])
    if feat:
        features.append(feat)
        targets.append(row["logS (mol/L)"])

# Train/test split
X_train, X_test, y_train, y_test = train_test_split(features, targets, test_size=0.2, random_state=42)

# Train model
model = RandomForestRegressor(n_estimators=100, random_state=42)
model.fit(X_train, y_train)

# Evaluate
preds = model.predict(X_test)
rmse = sqrt(mean_squared_error(y_test, preds))
print(f"Test RMSE: {rmse:.3f}")

# Save model
joblib.dump(model, "logS_model.pkl")