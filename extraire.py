"""Extrait les paramètres d'un ECG à partir de son image.

Usage : python extraire.py exemples/normal.png [autres images…]
Sans argument : génère puis analyse les 4 exemples fictifs.
"""
import sys

from visionecg import analyser, imprimer

EXEMPLES = {"normal": dict(fc=72), "tachycardie": dict(fc=115),
            "irregulier": dict(fc=95, irregularite=0.18), "qrs-large": dict(fc=64, largeur_qrs=2.4)}


def rapport(chemin):
    r = analyser(chemin)
    print(f"== {chemin}")
    print(f"   grille : {r['px_par_mm']:.2f} px/mm · étalonnage mesuré : {r['gain_mesure_mm_par_mv']:.1f} mm/mV")
    print(f"   {r['nb_battements']} battements · FC {r['fc']:.0f} /min ({r['frequence']}) · rythme {r['rythme']}")
    print(f"   RR moyen {r['rr_moyen_ms']:.0f} ms · SDNN {r['sdnn_ms']:.0f} ms")
    print(f"   onde R {r['amplitude_r_mv']:.2f} mV · QRS {r['qrs_ms']:.0f} ms ({r['qrs']})")


if __name__ == "__main__":
    fichiers = sys.argv[1:]
    if not fichiers:
        for nom, kw in EXEMPLES.items():
            imprimer(graine=2, **kw)[0].save(f"exemples/{nom}.png")
            fichiers.append(f"exemples/{nom}.png")
    for f in fichiers:
        rapport(f)
