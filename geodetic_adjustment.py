import numpy as np

# --- 1. JEODEZİK DÖNÜŞÜM FONKSİYONU ---
def geodetic_to_ecef(lat_deg, lat_min, lat_sec, lon_deg, lon_min, lon_sec, h):
    """Enlem, Boylam ve Elipsoid Yüksekliğini TUREF X,Y,Z'ye çevirir."""
    a = 6378137.0 
    f = 1 / 298.257223563 
    e2 = 2*f - f**2 
    lat = np.radians(lat_deg + lat_min/60.0 + lat_sec/3600.0)
    lon = np.radians(lon_deg + lon_min/60.0 + lon_sec/3600.0)
    N = a / np.sqrt(1 - e2 * np.sin(lat)**2)
    X = (N + h) * np.cos(lat) * np.cos(lon)
    Y = (N + h) * np.cos(lat) * np.sin(lon)
    Z = (N * (1 - e2) + h) * np.sin(lat)
    return np.array([X, Y, Z])

def perform_comprehensive_adjustment():
    # Terminal çıktı ayarları
    np.set_printoptions(precision=4, suppress=True, linewidth=200)

    # A. REFERANS NOKTALARI (TUREF PDF Verileri)
    ista_xyz = geodetic_to_ecef(41, 6, 16.01005, 29, 1, 9.62512, 147.2410)
    pala_xyz = geodetic_to_ecef(41, 5, 10.76464, 28, 57, 47.52407, 170.5481)
    
    fixed_coords = {'ISTA00TUR': ista_xyz, 'PALA': pala_xyz}

    # B. BAZ ÖLÇÜLERİ (6 Adet - B4 Hariç)
    baselines = [
        ('ISTA00TUR', '1701', -416.1682, 583.3060, 96.0653, [1.35792e-5, 8.3289e-6, 8.1942e-6, 6.0263e-6, 3.8325e-6, 7.6577e-6], 'B3'),
        ('ISTA00TUR', '1703', -258.5255, 862.3389, -248.3975, [1.25486e-5, 1.15957e-5, 2.29569e-5, 8.1907e-6, 1.04995e-5, 1.54703e-5], 'B1'),
        ('ISTA00TUR', '1702', -356.3439, 1036.5623, -273.0157, [5.2972e-6, 2.0774e-6, 2.2432e-6, 4.0283e-6, 1.7735e-6, 6.6627e-6], 'B2'),
        ('PALA', '1701', -3874.7995, 4059.1086, 1597.6165, [9.36615e-5, 5.71578e-5, 5.24039e-5, 4.05775e-5, 2.44522e-5, 5.13895e-5], 'B7'),
        ('PALA', '1703', -3717.1597, 4338.1358, 1253.1516, [9.48212e-5, 8.38954e-5, 1.507405e-4, 5.65865e-5, 7.37238e-5, 9.13276e-5], 'B5'),
        ('PALA', '1702', -3814.9799, 4512.3646, 1228.5268, [2.80600e-5, 1.09004e-5, 1.19072e-5, 2.05131e-5, 8.4779e-6, 3.33926e-5], 'B6')
    ]

    # Bilinmeyenler: 1701, 1702, 1703 (u = 9)
    unknown_points = ['1701', '1702', '1703']
    p_idx = {name: i for i, name in enumerate(unknown_points)}
    all_points = ['ISTA00TUR', 'PALA', '1701', '1702', '1703']

    # Yaklaşık Koordinatlar (Tablo hesabı için)
    approx_coords = {
        'ISTA00TUR': ista_xyz,
        'PALA': pala_xyz,
        '1701': ista_xyz + np.array([-416.1682, 583.3060, 96.0653]),
        '1702': ista_xyz + np.array([-356.3439, 1036.5623, -273.0157]),
        '1703': ista_xyz + np.array([-258.5255, 862.3389, -248.3975])
    }

    n_obs = len(baselines) * 3
    n_unkn = len(unknown_points) * 3
    A = np.zeros((n_obs, n_unkn))
    P = np.zeros((n_obs, n_obs))
    L = np.zeros(n_obs)

    # Matris Dolumu
    for i, (p_from, p_to, dx, dy, dz, cov, b_name) in enumerate(baselines):
        cov_m = np.array([[cov[0], cov[1], cov[3]], [cov[1], cov[2], cov[4]], [cov[3], cov[4], cov[5]]])
        P[i*3:i*3+3, i*3:i*3+3] = np.linalg.inv(cov_m)
        L[i*3:i*3+3] = [dx, dy, dz]
        if p_to in unknown_points: A[i*3:i*3+3, p_idx[p_to]*3:p_idx[p_to]*3+3] = np.eye(3)
        if p_from in unknown_points: A[i*3:i*3+3, p_idx[p_from]*3:p_idx[p_from]*3+3] = -np.eye(3)
        if p_from in fixed_coords: L[i*3:i*3+3] += fixed_coords[p_from]
        if p_to in fixed_coords: L[i*3:i*3+3] -= fixed_coords[p_to]

    # Çözüm: X = (A^T P A)^-1 * A^T P L
    N = A.T @ P @ A
    U = A.T @ P @ L
    Qx = np.linalg.inv(N)
    X_adj = Qx @ U
    
    # Rezidüeller ve İstatistikler
    V = A @ X_adj - L
    dof = n_obs - n_unkn  # 18 - 9 = 9
    m0_sq = (V.T @ P @ V) / dof
    Sigma_X = m0_sq * Qx

    # --- ÇIKTI RAPORU ---
    print("\n" + "█"*115)
    print("  TUREF DATUMUNDA DAYALI AĞ DENGELEME RAPORU (ISTA & PALA SABİT)")
    print("█"*115)

    print("\n[A] TASARIM MATRİSİ (A) [Boyut: 18x9]")
    print("-" * 115)
    # A matrisindeki 1, -1 ve 0'ların düzgün görünmesi için özel formatter
    np.set_printoptions(formatter={'float': lambda x: f"{x:5.0f}"})
    print(A)

    print("\n[L] GÖZLEMLER VE SABİT KOORDİNAT FARKLARI VEKTÖRÜ (L) [m]")
    print("-" * 115)
    np.set_printoptions(formatter={'float': lambda x: f"{x:15.4f}"})
    print(L.reshape(-1, 1))

    print("\n[1] DETAYLI AĞ İSTATİSTİKLERİ")
    print("-" * 55)
    print(f"Gözlem Sayısı (n)           : {n_obs}")
    print(f"Bilinmeyen Sayısı (u)       : {n_unkn}")
    print(f"Serbestlik Derecesi (f)     : {dof}")
    print(f"Varyans Faktörü (m0^2)      : {(m0_sq * 1e6):.6f} mm^2")
    print(f"Ağ Standart Sapması (m0)    : ±{np.sqrt(m0_sq)*1000:.2f} mm")

    print("\n[2] VARYANS-KOVARYANS MATRİSİ (Sigma_X) [mm^2]")
    print("-" * 115)
    np.set_printoptions(formatter={'float': lambda x: f"{x:10.2f}"})
    print(Sigma_X * 1e6)

    print("\n[3] BAZ ÖLÇÜLERİNE GETİRİLEN DÜZELTMELER (v = Ax - L) [mm]")
    print("-" * 75)
    print(f"{'BAZ':<6} | {'GÜZERGAH':<25} | {'v_X (mm)':>10} | {'v_Y (mm)':>10} | {'v_Z (mm)':>10}")
    print("-" * 75)
    for i, (p_from, p_to, dx, dy, dz, cov, b_name) in enumerate(baselines):
        v = V[i*3:i*3+3] * 1000
        print(f"{b_name:<6} | {p_from + ' -> ' + p_to:<25} | {v[0]:+10.2f} | {v[1]:+10.2f} | {v[2]:+10.2f}")

    print("\n[4] DENGELENMİŞ KESİN KOORDİNATLAR VE DOĞRULUKLAR")
    print("-" * 115)
    print(f"{'NOKTA':<12} | {'EKSEN':<5} | {'ÖNCE (m)':>15} | {'SONRA (m)':>15} | {'DÜZELTME (mm)':>15} | {'HATA (mm)':>10}")
    print("-" * 115)
    for pt in all_points:
        before = approx_coords[pt]
        if pt in fixed_coords:
            after = fixed_coords[pt]
            diff = np.zeros(3)
            err = np.zeros(3)
        else:
            after = X_adj[p_idx[pt]*3 : p_idx[pt]*3+3]
            diff = (after - before) * 1000
            err = np.sqrt(np.diag(Sigma_X)[p_idx[pt]*3 : p_idx[pt]*3+3]) * 1000

        axes = ['X', 'Y', 'Z']
        for i in range(3):
            label = pt if i == 0 else ""
            d_str = "0.00 (Sabit)" if pt in fixed_coords else f"{diff[i]:+15.2f}"
            e_str = "SABİT" if pt in fixed_coords else f"±{err[i]:.2f}"
            print(f"{label:<12} | {axes[i]:<5} | {before[i]:15.4f} | {after[i]:15.4f} | {d_str:>15} | {e_str:>10}")
        print("-" * 115)

    print("\n[P] AĞIRLIK MATRİSİ (P) [1/mm^2]")
    print("-" * 115)
    np.set_printoptions(formatter={'float': lambda x: f"{x:8.4f}" if abs(x) > 1e-10 else "  0.0000"})
    print(P / 1e6)
    print("█"*115)

if __name__ == "__main__":
    perform_comprehensive_adjustment()