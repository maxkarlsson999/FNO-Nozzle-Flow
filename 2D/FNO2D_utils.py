import torch
from torch.utils.data import DataLoader, TensorDataset
import numpy as np
import h5py



def data_loader(path_xy, path_TP, path_uv, batch_size):
    
    # ------------------------------------------------
    # 4) Data Loading from single HDF5 file
    # ------------------------------------------------

    with h5py.File(path_xy, 'r') as f:
        bcData = f['/BC'][:]    # shape [B, 7]
        X_all  = f['/X'][:].transpose(2,1,0)      # Shape [B, N_x, N_y]
        Y_all  = f['/Y'][:].transpose(2,1,0)      # Shape [B, N_x, N_y]

    # --- load T, P  ---
    with h5py.File(path_TP, 'r') as f:
        T_all = f['/T'][:].transpose(2,1,0)      # Shape [B, N_x, N_y]
        P_all = f['/P'][:].transpose(2,1,0)      # Shape [B, N_x, N_y]

    # --- load u, v  ---
    with h5py.File(path_uv, 'r') as f:
        u_all = f['/u'][:].transpose(2,1,0)      # Shape [B, N_x, N_y]
        v_all = f['/v'][:].transpose(2,1,0)      # Shape [B, N_x, N_y]
        vel_mag = np.sqrt(u_all**2 + v_all**2)       # Convert to velocity magnitude

    # ------------------------------------------------
    # Preallocate arrays
    # ------------------------------------------------
    N = bcData.shape[0]
    N, N_x, N_y = X_all.shape 


    # 1) Inputs: T0_L, P0_L, P_R  -> shape [B, 3]
    inputs = bcData[:, 0:3].astype(np.float32)  # pick out the first 3 columns

    # 2) Outputs: (T, P, u) each [B, N_x, N_y], combine into [B, N_x, N_y, 3]
    outputs = np.empty((N, N_x, N_y, 3), dtype=np.float32)
    outputs[:, :, :, 0] = T_all.astype(np.float32)
    outputs[:, :, :, 1] = P_all.astype(np.float32)
    outputs[:, :, :, 2] = vel_mag.astype(np.float32)

    # 3) x_grid, y_grid each [B, N_x, N_y]
    x_grid = X_all.astype(np.float32)
    y_grid = Y_all.astype(np.float32)

    # ------------------------------------------------
    # Normalize inputs T0, P0, P_B
    # ------------------------------------------------
    inp_max = np.max(inputs, axis=0)
    inp_min = np.min(inputs, axis=0)
    inputs_norm = (inputs - inp_min) / ( (inp_max - inp_min) )

    # ------------------------------------------------
    # Normalize outputs (T, P, u)
    # ------------------------------------------------
    out_max = np.max(outputs, axis=(0,1,2))
    out_min = np.min(outputs, axis=(0,1,2))
    outputs_norm = (outputs - out_min) / ( (out_max - out_min) )

    # ------------------------------------------------
    # Normalize the x and y grids
    # ------------------------------------------------
    x_max = np.max(x_grid, axis=(0,1,2))
    x_min = np.min(x_grid, axis=(0,1,2))
    x_norm = (x_grid - x_min) / ( (x_max - x_min) )

    y_max = np.max(y_grid, axis=(0,1,2))
    y_min = np.min(y_grid, axis=(0,1,2))
    y_norm = (y_grid - y_min) / ( (y_max - y_min)  )

    # ------------------------------------------------
    # Convert to torch tensors, create Dataset
    # ------------------------------------------------
    inputs_torch  = torch.tensor(inputs_norm,  dtype=torch.float32)
    outputs_torch = torch.tensor(outputs_norm, dtype=torch.float32)
    x_torch       = torch.tensor(x_norm,       dtype=torch.float32)
    y_torch       = torch.tensor(y_norm,       dtype=torch.float32)

    # Remove unphysical samples before loading
    mask = torch.ones(N, dtype=torch.bool)
    idx_remove = np.array([3295, 1537, 3803, 1299, 3877, 1646,  961, 2832,  602, 4444, 4494, 2321])
    mask[idx_remove] = False

    # 2) Index each tensor with that mask
    new_inputs  = inputs_torch[mask]    
    new_outputs = outputs_torch[mask]
    new_x       = x_torch[mask]
    new_y       = y_torch[mask]

    dataset = TensorDataset(new_inputs, new_outputs, new_x, new_y)

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
    Prediction tensor, shape: (Batch, N_x, N_y, Output Channels)

    Returns:
    Target tensor, shape: (Output Channels)
    """
    split_loss =  (torch.mean((pred-target)**2, dim=[0,1,2])) /(torch.mean(target**2, dim=[0,1,2]))
    return torch.sqrt(split_loss)

def relative_l2(pred, target):
    """
    Computes relative RMSE error normalized by targets, average taken over each case and grid

    Args: 
    Prediction tensor, shape: (Batch, N_x, N_y, Output Channels)

    Returns:
    Target tensor, shape: (Output Channels)
    """

    split_loss = torch.mean(((pred-target)**2) /(target**2 +0.01), dim=[0,1,2])
    return torch.sqrt(split_loss)


def relative_min_max(pred, target):
    """
    Computes relative RMSE error normalized by max-min of target, average taken over each case and grid

    Args: 
    Prediction tensor, shape: (Batch, N_x, N_y, Output Channels)

    Returns:
    Target tensor, shape: (Output Channels)
    """

    split_loss = torch.mean(((pred-target)**2) /((torch.amax(target, dim=[1,2], keepdim=True) - torch.amin(target, dim=[1,2], keepdim=True))**2), dim=[0,1,2])
    return torch.sqrt(split_loss)

def normal_RMSE(pred, target): 
    """
    Computes regular RMSE over batch and grid

    Args: 
    Prediction tensor, shape: (Batch, N_x, N_y, Output Channels)

    Returns:
    Target tensor, shape: (Output Channels)
    """
    split_loss = torch.mean((pred-target)**2, dim=[0,1,2])
    return torch.sqrt(split_loss)
