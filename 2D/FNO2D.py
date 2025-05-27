import torch
import torch.nn as nn
import torch.nn.functional as F

class SpectralConv2d(nn.Module):
    """
    A 2D spectral convolution using truncated Fourier transforms.
    
    Steps:
        1. rFFT in 2D
        2. Multiply low-frequency modes by learnable complex weights
        3. iFFT to get back to real space
    """
    
    def __init__(self, in_channels, out_channels, modes1, modes2):
        super().__init__()
        self.in_channels = in_channels
        self.out_channels = out_channels
        self.modes1 = modes1  
        self.modes2 = modes2  
    
        scale = 1.0 / (in_channels * out_channels)
        self.weights = nn.Parameter(
            scale * torch.randn(in_channels, out_channels, modes1, modes2, dtype=torch.cfloat)
        )
    
    def _compl_matmul_2d(self, x, w):
        """
        Forward pass of spectral convolution
    
        Args:
        x: Input tensor, shape:(batch, in_channels, N_x, N_y)
        
        Returns: 
        x_out: Output tensor, shape: (batch, out_channels, N_x, N_y) in real space.
        """
        return torch.einsum("bixy,ioxy->boxy", x, w)
    
    def forward(self, x: torch.Tensor):
        """
        Forward pass of spectral convolution
    
        Args:
        x: Input tensor, shape:(batch, in_channels, N_x, N_y)
        
        Returns: 
        x_out: Output tensor, shape: (batch, out_channels, N_x, N_y) in real space.
        """
        B, C, N_x, N_y = x.shape
    
        x_ft = torch.fft.rfft2(x, norm="ortho")  
    
        out_ft = torch.zeros(
            B, self.out_channels, N_x, N_y // 2 + 1, dtype=torch.cfloat, device=x.device
        )
        out_ft[:, :, : self.modes1, : self.modes2] = self._compl_matmul_2d(
            x_ft[:, :, : self.modes1, : self.modes2], self.weights
        )
    
        x_out = torch.fft.irfft2(out_ft, s=(N_x, N_y), norm="ortho")  # (B, out_channels, N_x, N_y)
        return x_out


class FNO2d(nn.Module):
    """
    2D Fourier Neural Operator model.

    Lifting layer -> Spectral convolution layers + Pointwise convolutions -> Projection layers
    """

    def __init__(self, modes1, modes2, width, hidden_mlp, N_x, N_y, N_fourier_layers):
        super().__init__()
        self.modes1 = modes1
        self.modes2 = modes2
        self.width = width
        self.N_x = N_x  
        self.N_y = N_y  
        self.N_fourier_layers = N_fourier_layers

        # --- Lifting Layer ---
        self.mlp1 = nn.Linear(5, width)

        # --- Spectral Conv Layers + Pointwise Conv ---
        self.conv_layers = nn.ModuleList(
            [SpectralConv2d(width, width, modes1, modes2) for _ in range(N_fourier_layers)]
        )
        self.w_layers = nn.ModuleList(
            [nn.Conv2d(width, width, kernel_size=1) for _ in range(N_fourier_layers)]
        )

        # --- Multi-Head Decoders ---
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
        Forward pass of FNO2D 

        Args: 
        bc: Initial conditions tensor, shape: (batch, 3)
        x: x-coordinates, shape: (batch, N_x, N_y)
        y: y-coordinates, shape: (batch, N_x, N_y)

        Returns: 
        final_out: Output fields of temperature, pressure and velocity, shape: (batch, N_x, N_y, 3)
        """
        B, N_x, N_y = x.shape

        bc_exp = bc.view(B, 1, 1, 3).expand(-1, N_x, N_y, -1)  
        xy = torch.stack([x, y], dim=-1)                   
        inp = torch.cat([bc_exp, xy], dim=-1)              

        h = self.mlp1(inp)                                
        h = h.permute(0, 3, 1, 2)                        

        # ------------- spectral blocks -------------
        for conv, w in zip(self.conv_layers, self.w_layers):
            h = h + F.gelu(conv(h) + w(h))                

        h = h.permute(0, 2, 3, 1).contiguous()           
        h_flat = h.view(B * N_x * N_y, self.width)        

        # Apply decoders separately
        T = self.temp_decoder(h_flat).view(B, N_x, N_y, 1)
        P = self.pressure_decoder(h_flat).view(B, N_x, N_y, 1)
        u = self.vel_decoder(h_flat).view(B, N_x, N_y, 1)
            # Concatenate outputs back into single tensor: [T, P, u]
        return torch.cat([T, P, u], dim=-1)               # (B, N_x, N_x, 3)
