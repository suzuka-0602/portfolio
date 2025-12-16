# portfolio
自己紹介のportfolioを作成しました．

## MATLAB でのローレンツフィッティング

`matlab/lorentzian_fitting.m` と `matlab/lorentzian_regions.m` を使って、各 `.spe`
スペクトルに対するローレンツフィッティングを行い、ピーク強度 `A_1` などの係数を
構造体にまとめて `lorentzian_results.mat` に保存します。`matlab/Caliblation.m`
でスペクトル行番号と cm^-1 の対応付けを設定してください（`wn_calib = Caliblation(x)` の
形で利用できます）。

1. `matlab/lorentzian_regions.m` にフィッティング対象のスペクトル名、中心
   wavenumber、ウィンドウ幅、ピーク数（`peak_centers_cm1`）を記述します。
2. `matlab/Caliblation.m` で wavenumber 軸を定義します（線形間隔・多項式など自由に記述）。
3. MATLAB で `lorentzian_fitting('スペクトルフォルダ')` を実行すると、各スペクトルの
   `A_k` が `results.<スペクトル名>.lorentz_A_k` として保存されます。
