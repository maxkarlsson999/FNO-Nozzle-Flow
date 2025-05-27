import torch
import torch.nn as nn
import torch.nn.functional as F

class SpectralConv1d(nn.Module):
    """
    A 1D spectral convolution using truncated Fourier transforms.

    Steps:
      1. rFFT in 1D
      2. Multiply low-frequency modes by learnable complex weights
      3. iFFT to get back to real space
    """
    def __init__(self, in_channels, out_channels, modes):
        """
        Args:
        in_channels : # of input channels
        out_channels: # of output channels
        modes       : # of low-frequency modes to keep
        """
        super(SpectralConv1d, self).__init__()
        self.in_channels = in_channels
        self.out_channels = out_channels
        self.modes = modes
        self.scale = 1.0 / (in_channels * out_channels)
        self.weights = nn.Parameter(
            self.scale * torch.randn(in_channels, out_channels, modes, dtype=torch.cfloat)
        )

    def complex_matmul_1d(self, x, w):
        """
        Complex tensor multiplication
        
        Args:
        x: (batch, in_channels, modes) 
        w: (in_channels, out_channels, modes)

        Returns:
        Tensor of shape (batch, out_channels, modes)
        """
        return torch.einsum("bik,iok->bok", x, w)

    def forward(self, x):
        """
        Forward pass of spectral convolution

        Args:
        x: Input tensor, shape:(batch, in_channels, N_x)
        
        Returns: 
        x_out: Output tensor, shape: (batch, out_channels, N_x) in real space.
        """
        B, C, N_x = x.shape
        x_ft = torch.fft.rfft(x, norm="ortho")  # (B, C, N_x//2+1)
        out_ft = torch.zeros(
            (B, self.out_channels, N_x//2 + 1), dtype=torch.cfloat, device=x.device
        )
        freq_slice = slice(0, self.modes)
        out_ft[:, :, freq_slice] = self.complex_matmul_1d(x_ft[:, :, freq_slice], self.weights)
        x_out = torch.fft.irfft(out_ft, n=N_x, norm="ortho")
        return x_out

class FNO1d(nn.Module):
    """
    1D Fourier Neural Operator model.

    Lifting layer -> Spectral convolution layers + Pointwise convolutions -> Projection layers
    """
    
    def __init__(self, modes, width, hidden_mlp, N_x, N_fourier_layers):
        super(FNO1d, self).__init__()
        self.modes = modes
        self.width = width
        self.N_x = N_x
        self.N_fourier_layers = N_fourier_layers

        # --- Lifting Layer ---
        self.mlp1 = nn.Linear(5, width)

        # --- Spectral Conv Layers + Pointwise Conv ---
        self.conv_layers = nn.ModuleList(
            [SpectralConv1d(width, width, modes) for _ in range(N_fourier_layers)]
        )
        self.w_layers = nn.ModuleList(
            [nn.Conv1d(width, width, kernel_size=1) for _ in range(N_fourier_layers)]
        )

        # --- Multi-Head Projection Layers ---
        # Temperature decoder
        self.temp_decoder = nn.Sequential(
            nn.Linear(width, width),
            nn.GELU(),
            nn.Linear(width, hidden_mlp),
            nn.GELU(),
            nn.Linear(hidden_mlp, 1) 
        )

        # Pressure decoder
        self.pressure_decoder = nn.Sequential(
            nn.Linear(width, width),
            nn.GELU(),
            nn.Linear(width, hidden_mlp),
            nn.GELU(),
            nn.Linear(hidden_mlp, 1)  
        )
        
        # Velocity decoder
        self.vel_decoder = nn.Sequential(
            nn.Linear(width, width),
            nn.GELU(),
            nn.Linear(width, hidden_mlp),
            nn.GELU(),
            nn.Linear(hidden_mlp, 1)  
        )

    def forward(self, bc, x, y):
        """
        Forward pass of FNO1D 

        Args: 
        bc: Initial conditions tensor, shape: (batch, 3)
        x: x-coordinates, shape: (batch, N_x)
        y: Cross-sectional area, shape: (batch, N_x)

        Returns: 
        final_out: Output fields of temperature, pressure and velocity, shape: (batch, N_x, 3)
        """
        
        bc_expanded = bc.unsqueeze(1).repeat(1, self.N_x, 1)
        xy = torch.stack([x, y], dim=-1)
        inp = torch.cat([bc_expanded, xy], dim=-1)

        x = self.mlp1(inp)  # (B, N_x, width)
        x = x.permute(0, 2, 1)  

        for conv, w in zip(self.conv_layers, self.w_layers):
            x = x + F.gelu(conv(x) + w(x))
        x = x.permute(0, 2, 1)  # (B, N_x, width)


        # Apply decoders separately
        temp_out = self.temp_decoder(x)           # (B, N_x, 1)
        pressure_out = self.pressure_decoder(x) # (B, N_x, 1)
        vel_out = self.vel_decoder(x)           # (B, N_x, 1)
        # Concatenate outputs back into single tensor: [T, P, u]
        final_out = torch.cat([temp_out, pressure_out, vel_out], dim=-1)  # (B, N_x, 3)

        return final_out
