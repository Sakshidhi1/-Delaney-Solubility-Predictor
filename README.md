#🧪  Delaney-Solubility-Predictor
This project is a *web application* that predicts the aqueous solubility of molecules (logS in mol/L) using a *Random Forest model* trained on the Delaney dataset. It uses *SMILES* as input and outputs a predicted solubility value.

## 🎥 Demo Video
Watch the app in action:  

https://github.com/user-attachments/assets/c226eaf8-d79a-4ffc-98d8-f75a4ac38749

## 🚀 Features

- 🔬 Predict aqueous solubility from a SMILES string
- 🧠 Uses a trained machine learning model (Random Forest)
- 🧾 Based on *Delaney 2004 dataset*
- 🌐 Built using *Streamlit, **scikit-learn, **RDKit*

  # 📁 Project Structure

delaney_solubility_app/ ├── delaney.csv            # Dataset file ├── train_model.py         # Script to train the model ├── app.py                 # Streamlit app ├── logS_model.pkl         # Trained model (generated) ├── README.md              # This file └── venv/                  # Virtual environment (optional)

## ⚙️ Installation Instructions

### 1️⃣ Clone or Download

Create a folder delaney_solubility_app and place all files inside.

### 2️⃣ Create Virtual Environment

```bash
cd Desktop\delaney_solubility_app
python -m venv venv
venv\Scripts\activate


3️⃣ Install Dependencies

pip install pandas scikit-learn streamlit joblib rdkit-pypi


🏋️‍♂️ Train the Model

python train_model.py

This will generate logS_model.pkl and print the test RMSE.



🚀 Run the Web App

streamlit run app.py

Then open your browser at http://localhost:8501



🧪 Example SMILES Inputs
Ethanol	CCO
Benzene	C1=CC=CC=C1
Aspirin	CC(=O)OC1=CC=CC=C1C(=O)O
Caffeine	CN1C=NC2=C1C(=O)N(C(=O)N2C)C


## 📊 Dataset Source

Delaney, J. S. (2004). ESOL: Estimating aqueous solubility directly from molecular structure.
J. Chem. Inf. Comput. Sci., 44(3), 1000–1005
Dataset: DeepChem GitHub – delaney-processed.csv
Dataset kagggle- delaney.csv


## ✅ Notes

Your dataset must include these columns:

SMILES

logS (mol/L)


If your column names differ, update them in train_model.py.

📸 Screenshots (Optional)
![WhatsApp Image 2025-07-09 at 7 40 26 PM](https://github.com/user-attachments/assets/7b739ea2-bfe1-4f6b-878c-c33735a73a6d)






