# lstm_sysid.py
import numpy as np
import pandas as pd
import tensorflow as tf
from tensorflow.keras.models import Sequential, load_model
from tensorflow.keras.layers import LSTM, Dense
from tensorflow.keras.callbacks import EarlyStopping

# ---------------------- Helper Functions ---------------------- #
def createInputData(u, th, th0=0):
    shiftedOutput = [0] * (len(u) + 1)
    shiftedOutput[0] = th0
    shiftedOutput[1:] = th[:-1]
    combinedInput = np.transpose([u, shiftedOutput])
    return np.array(combinedInput), np.array(shiftedOutput[:-1])

def create_sequences(data, labels, seq_length):
    X, y = [], []
    for i in range(len(data) - seq_length):
        # print(i)
        X.append(data[i:i + seq_length])
        y.append(labels[i + seq_length])
    return np.array(X), np.array(y)

# ---------------------- Main Functions ---------------------- #
def train_model(file_path, save_path, seq_length=10, epochs=100):

    # Load data
    data = pd.read_csv(file_path)
    u = data.iloc[:,0].values
    th = data.iloc[:,1].values

    # Create input sequences
    combinedInput, shiftedOutput = createInputData(u, th)
    X, y = create_sequences(combinedInput, th, seq_length)
    X = X.reshape((X.shape[0], X.shape[1], X.shape[2]))

    # Normalize
    mean_X, std_X = X.mean(axis=(0,1), keepdims=True), X.std(axis=(0,1), keepdims=True)
    X = (X - mean_X) / std_X
    mean_y, std_y = y.mean(), y.std()
    y = (y - mean_y) / std_y

    # Chronological train/val split
    split_idx = int(0.8 * len(X))
    X_train, X_val = X[:split_idx], X[split_idx:]
    y_train, y_val = y[:split_idx], y[split_idx:]

    # Define model
    model = Sequential()
    model.add(LSTM(64, activation='tanh', input_shape=(seq_length, X.shape[2]), return_sequences=True))
    model.add(LSTM(64, activation='tanh'))
    model.add(Dense(1))

    model.compile(optimizer=tf.keras.optimizers.Adam(learning_rate=0.0005, clipnorm=1.0),
                  loss='mse')

    # Train
    early_stopping = EarlyStopping(monitor='val_loss', patience=10, restore_best_weights=True)
    history = model.fit(X_train, y_train,
                        validation_data=(X_val, y_val),
                        epochs=epochs,
                        batch_size=16,
                        callbacks=[early_stopping],
                        verbose=0)

    # Predict on validation
    y_pred = model.predict(X_val)
    y_pred_denorm = y_pred.flatten() * std_y + mean_y
    y_val_denorm = y_val.flatten() * std_y + mean_y
    rmse = np.sqrt(np.mean((y_pred_denorm - y_val_denorm)**2))

    # Save model and normalization parameters
    model.save(save_path)
    np.savez(save_path + "_norm.npz", mean_X=mean_X, std_X=std_X, mean_y=mean_y, std_y=std_y)

    return float(rmse)

def predict_model(model_path, input_array):

    # Load model and normalization
    model = load_model(model_path)
    norm = np.load(model_path + "_norm.npz", allow_pickle=True)
    mean_X, std_X = norm['mean_X'], norm['std_X']
    mean_y, std_y = float(norm['mean_y']), float(norm['std_y'])

    X = np.array(input_array, dtype=float)
    X = (X - mean_X) / std_X

    y_pred = model.predict(X)
    y_pred_denorm = y_pred.flatten() * std_y + mean_y
    return y_pred_denorm.tolist()


# train_model('file.csv', 'lstm_model_internship.keras', 10, 50)
