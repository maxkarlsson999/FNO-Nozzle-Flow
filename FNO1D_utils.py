import torch
from torch.utils.data import DataLoader, TensorDataset
import numpy as np
import h5py

def data_loader(path, batch_size):
    
    # ------------------------------------------------
    # 4) Data Loading from single HDF5 file
    # ------------------------------------------------
    
    # Open HDF5 file 
    with h5py.File(path, 'r') as f:
        # Read BC data, shape = [B, 7] 
        bcData = f['BC'][:]  # columns: [T0, P0, P_B, x_throat, r_inlet, r_throat, r_out]
        # Read geometry and solution arrays, each shape = [B, N_x]
        X_all = f['X'][:]    # x-coordinates
        Y_all = f['Y'][:]    # nozzle radius
        Y_all = (Y_all**2)*np.pi # Convert to cross-sectional area
        T_all = f['T'][:]    # temperature
        P_all = f['P'][:]    # pressure
        u_all = f['u'][:]    # velocity

    # ------------------------------------------------
    # Preallocate arrays
    # ------------------------------------------------
    B = bcData.shape[0]
    N_x = X_all.shape[1]   

    # 1) Inputs: T0, P0, P_B  -> shape [B, 3]
    inputs = bcData[:, 0:3].astype(np.float32)  # pick out the first 3 columns

    # 2) Outputs: (T, P, u) each [B, N_x], combine into [B, N_x, 3]
    outputs = np.empty((B, N_x, 3), dtype=np.float32)
    outputs[:, :, 0] = T_all.astype(np.float32)
    outputs[:, :, 1] = P_all.astype(np.float32)
    outputs[:, :, 2] = u_all.astype(np.float32)

    # 3) x_grid, y each [B, N_x]
    x_grid = X_all.astype(np.float32)
    y_grid = Y_all.astype(np.float32)

    # ------------------------------------------------
    # Normalize inputs T0, P0, P_B
    # ------------------------------------------------
    inp_max = np.max(inputs, axis=0)
    inp_min = np.min(inputs, axis=0)
    inputs_norm = (inputs - inp_min) / ( (inp_max - inp_min)) 

    # ------------------------------------------------
    # Normalize outputs (T, P, u)
    # ------------------------------------------------
    out_max = np.max(outputs, axis=(0,1))
    out_min = np.min(outputs, axis=(0,1))
    outputs_norm = (outputs - out_min) / ( (out_max - out_min) ) 

    # ------------------------------------------------
    # Normalize x and y
    # ------------------------------------------------
    x_max = np.max(x_grid, axis=(0,1))
    x_min = np.min(x_grid, axis=(0,1))
    x_norm = (x_grid - x_min) / ((x_max - x_min))

    y_max = np.max(y_grid, axis=(0,1))
    y_min = np.min(y_grid, axis=(0,1))
    y_norm = (y_grid - y_min) / ((y_max - y_min)) 

    # ------------------------------------------------
    # Convert to torch tensors, create Dataset
    # ------------------------------------------------
    inputs_torch  = torch.tensor(inputs_norm,  dtype=torch.float32)
    outputs_torch = torch.tensor(outputs_norm, dtype=torch.float32)
    x_torch       = torch.tensor(x_norm,       dtype=torch.float32)
    y_torch       = torch.tensor(y_norm,       dtype=torch.float32)

    dataset = TensorDataset(inputs_torch, outputs_torch, x_torch, y_torch)

    # ------------------------------------------------
    # Train/test split (80/20), create DataLoaders
    # ------------------------------------------------
    train_size = int(0.8 * len(dataset))
    test_size  = len(dataset) - train_size
    train_dataset, test_dataset = torch.utils.data.random_split(dataset, [train_size, test_size])

    train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
    test_loader  = DataLoader(test_dataset,  batch_size=batch_size, shuffle=False)
    return train_loader, test_loader, inp_max, inp_min, out_max, out_min, x_max, x_min, y_max, y_min


def relative_l2_separate_mean(pred, target):
    """
    Computes relative RMSE error normalized by average of targets

    Args: 
    Prediction tensor, shape: (Batch, Grid, Output Channels)

    Returns:
    Target tensor, shape: (Output Channels)
    """
    
    split_loss =  (torch.mean((pred-target)**2, dim=[0,1])) /(torch.mean(target**2, dim=[0,1]))
    return torch.sqrt(split_loss)

def relative_l2(pred, target):
    """
    Computes relative RMSE error normalized by targets, average taken over each case and grid
    
    Args: 
    Prediction tensor, shape: (Batch, Grid, Output Channels)

    Returns:
    Target tensor, shape: (Output Channels)
    """
    
    split_loss = torch.mean(((pred-target)**2) /(target**2 +0.01), dim=[0,1]) # Prevent division by 0 using 0.01
    return torch.sqrt(split_loss)


def relative_min_max(pred, target):
    """
    Computes relative RMSE error normalized by max-min of target, average taken over each case and grid
    
    Args: 
    Prediction tensor, shape: (Batch, Grid, Output Channels)

    Returns:
    Target tensor, shape: (Output Channels)
    """

    split_loss = torch.mean(((pred-target)**2) /((torch.amax(target, dim=[1], keepdim=True) - torch.amin(target, dim=[1], keepdim=True))**2), dim=[0,1])
    return torch.sqrt(split_loss)

def normal_RMSE(pred, target): 
    """
    Computes regular RMSE over batch and grid
    
    Args: 
    Prediction tensor, shape: (Batch, Grid, Output Channels)

    Returns:
    Target tensor, shape: (Output Channels)
    """
    
    split_loss = torch.mean((pred-target)**2, dim=[0,1])
    return torch.sqrt(split_loss)

def relative_min_max_MSE(pred, target):
    """
    Computes relative MSE error normalized by max-min of target, average taken over each case and grid
    
    Args: 
    Prediction tensor, shape: (Batch, Grid, Output Channels)

    Returns:
    Target tensor, shape: (Output Channels)
    """

    split_loss = torch.mean(((pred-target)**2) /((torch.amax(target, dim=[1], keepdim=True) - torch.amin(target, dim=[1], keepdim=True))**2), dim=[0,1])
    return split_loss
    