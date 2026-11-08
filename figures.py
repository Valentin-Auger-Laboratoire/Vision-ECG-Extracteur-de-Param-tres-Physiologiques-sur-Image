"""Génère images/vision-ecg.png (nécessite matplotlib) : python figures.py"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

from visionecg import imprimer, analyser

im, v = imprimer(fc=95, irregularite=0.18, graine=2, duree=6)
r = analyser(im)
tt, vv = r["signal"]
fig, ax = plt.subplots(2, 1, figsize=(12, 5.4), gridspec_kw={"height_ratios": [1.1, 1]})
ax[0].imshow(np.asarray(im)); ax[0].axis("off")
ax[0].set_title("Image d'entrée : ECG imprimé sur papier millimétré (fictif)")
ax[1].plot(tt, vv, "k", lw=1, label="tracé numérisé depuis l'image")
ax[1].plot(r["instants_r"], np.interp(r["instants_r"], tt, vv), "o", c="#c8102e", label="pics R détectés")
for a, b in zip(r["instants_r"][:-1], r["instants_r"][1:]):
    ax[1].annotate(f"{(b - a) * 1000:.0f}", ((a + b) / 2, 1.35), ha="center", fontsize=7, color="#555")
ax[1].set(xlabel="temps (s)", ylabel="mV", ylim=(-0.5, 1.6),
          title=f"FC {r['fc']:.0f} /min · rythme {r['rythme']} · SDNN {r['sdnn_ms']:.0f} ms · R {r['amplitude_r_mv']:.2f} mV · QRS {r['qrs_ms']:.0f} ms  (RR en ms)")
ax[1].grid(alpha=.3); ax[1].legend(loc="lower right", fontsize=8)
fig.tight_layout(); fig.savefig("images/vision-ecg.png", dpi=120)
