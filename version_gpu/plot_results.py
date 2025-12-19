#!/usr/bin/env python3
"""
Visualize Rd and Tt data from MCML GPU simulations
"""
import pandas as pd
import matplotlib.pyplot as plt
import os
import re

def load_and_plot(csv_filename, thickness, color, marker):
    filepath =  f'output_gpu/{csv_filename}'
    print(f"Processing {filepath}...")
    
    # Pre-clean file in memory
    try:
        data_lines = []
        if os.path.exists(filepath):
            with open(filepath, 'r') as f:
                for line in f:
                    if 'mco' in line or 'output,Rd,Tt' in line:
                        data_lines.append(line)
        
        if not data_lines:
            print(f"Warning: {filepath} is empty or invalid.")
            return None

        # Write to temp file
        with open('temp_clean.csv', 'w') as f:
            f.writelines(data_lines)
            
        df = pd.read_csv('temp_clean.csv')
        
        # Extract wavelength
        if 'wavelength' not in df.columns:
             df['wavelength'] = df['output'].apply(lambda x: float(re.findall(r'\d+', str(x))[0]) if re.findall(r'\d+', str(x)) else float('nan'))
        
        df = df.sort_values('wavelength')
        
        print(f"  {thickness}: {len(df)} wavelengths: {df['wavelength'].min():.0f} - {df['wavelength'].max():.0f} nm")
        
        # Plot Rd
        plt.figure(figsize=(10, 6))
        plt.plot(df['wavelength'], df['Rd'], color=color, linewidth=1.0, marker=marker, markersize=2, label=f'Rd ({thickness})')
        plt.title(f'Diffuse Reflectance (Rd) - {thickness}', fontsize=16)
        plt.xlabel('Wavelength (nm)')
        plt.ylabel('Rd')
        plt.grid(True, alpha=0.3)
        plt.legend()
        plt.savefig(f'output_gpu/Rd_{thickness}.png', dpi=300)
        plt.close()
        
        # Plot Tt
        plt.figure(figsize=(10, 6))
        plt.plot(df['wavelength'], df['Tt'], color=color, linewidth=1.0, marker=marker, markersize=2, label=f'Tt ({thickness})')
        plt.title(f'Transmittance (Tt) - {thickness}', fontsize=16)
        plt.xlabel('Wavelength (nm)')
        plt.ylabel('Tt')
        plt.grid(True, alpha=0.3)
        plt.legend()
        plt.savefig(f'output_gpu/Tt_{thickness}.png', dpi=300)
        plt.close()
        
        return df
        
    except Exception as e:
        print(f"Error: {e}")
        return None

def plot_combined():
    configs = [('summary_3mm.csv', '3mm', 'blue', 'o'),
               ('summary_4mm.csv', '4mm', 'green', 's'),
               ('summary_5mm.csv', '5mm', 'red', '^')]
    
    # Store dataframes
    dfs = {}
    
    for csv, thick, color, marker in configs:
        filepath = f'output_gpu/{csv}'
        try:
            # Pre-clean file in memory
            data_lines = []
            if os.path.exists(filepath):
                with open(filepath, 'r') as f:
                    for line in f:
                        if 'mco' in line or 'output,Rd,Tt' in line:
                            data_lines.append(line)
                            
            if not data_lines:
                continue

            # Parse in memory
            from io import StringIO
            df = pd.read_csv(StringIO(''.join(data_lines)))
            
            # Extract wavelength
            if 'wavelength' not in df.columns:
                    df['wavelength'] = df['output'].apply(lambda x: float(re.findall(r'\d+', str(x))[0]) if re.findall(r'\d+', str(x)) else float('nan'))
            
            df = df.sort_values('wavelength')
            dfs[thick] = (df, color, marker)
        except Exception as e:
            print(f"Error reading {csv}: {e}")

    # Plot Combined Rd
    plt.figure(figsize=(10, 6))
    for thick, (df, color, marker) in dfs.items():
        plt.plot(df['wavelength'], df['Rd'], color=color, linewidth=1.0, marker=marker, markersize=2, label=f'Rd ({thick})')
    
    plt.title('Diffuse Reflectance (Rd) - Combined', fontsize=16)
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Rd')
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig('output_gpu/Rd_combined.png', dpi=300)
    print("Saved output_gpu/Rd_combined.png")
    plt.close()
    
    # Plot Combined Tt
    plt.figure(figsize=(10, 6))
    for thick, (df, color, marker) in dfs.items():
        plt.plot(df['wavelength'], df['Tt'], color=color, linewidth=1.0, marker=marker, markersize=2, label=f'Tt ({thick})')
    
    plt.title('Transmittance (Tt) - Combined', fontsize=16)
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Tt')
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig('output_gpu/Tt_combined.png', dpi=300)
    print("Saved output_gpu/Tt_combined.png")
    plt.close()

if __name__ == "__main__":
    configs = [('summary_3mm.csv', '3mm', 'blue', 'o'),
               ('summary_4mm.csv', '4mm', 'green', 's'),
               ('summary_5mm.csv', '5mm', 'red', '^')]
               
    for csv, thick, color, marker in configs:
        load_and_plot(csv, thick, color, marker)
        
    if os.path.exists('temp_clean.csv'):
        os.remove('temp_clean.csv')
    
    print("-" * 30)
    print("Generating combined plots...")
    plot_combined()
    
    print("Done.")
