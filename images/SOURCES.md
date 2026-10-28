# Image sources

Every image here is our own work, in the public domain, or under a Creative Commons licence that allows reuse with attribution.
Each file is at most 800 px wide and 300 KB.
The photos were downloaded from Wikimedia Commons as the 960 px rendition, resized to 800 px wide (Lanczos) and re-encoded as JPEG quality 85; nothing else was changed.

| File | Purpose | Author, licence, source |
|---|---|---|
| `synthetic.png` | Easiest case: smooth sky and ground, two solid objects, fully deterministic | Drawn for this project (200x120): a vertical sky gradient, green ground below row 90, a red sun (centre (50, 40), radius 18) and a dark grey house (x 120-159, y 50-89) with a brown door. Carving to 50% width should keep the sun and the house intact. |
| `acacia.jpg` | Best case: one compact subject on a wide, plain sky and plain dry grass | Francisco Anzola, *Lone Tree*. [CC BY 3.0](https://creativecommons.org/licenses/by/3.0/). [Wikimedia Commons](https://commons.wikimedia.org/wiki/File:Lone_Tree_(250844427).jpeg). |
| `sailboat.jpg` | Best case: a thin, sharp silhouette on a smooth sunset sky and calm sea | Dietmar Rabich, *Zoutelande (NL), Strand, Blick auf die Nordsee* (2022). [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/). [Wikimedia Commons](https://commons.wikimedia.org/wiki/File:Zoutelande_(NL),_Strand,_Blick_auf_die_Nordsee_--_2022_--_4984.jpg). The resized copy is shared under the same licence. |
| `broadway_tower.jpg` | The classic seam-carving example (Avidan and Shamir 2007, Wikipedia, the assignment brief): a building and a small person on sky and grass | Newton2, cropped by Yummifruitbat, *Broadway tower edit*. [CC BY 2.5](https://creativecommons.org/licenses/by/2.5/). [Wikimedia Commons](https://commons.wikimedia.org/wiki/File:Broadway_tower_edit.jpg). |
| `starry_night.jpg` | Hardest case: busy texture everywhere, so every seam has to cut through detail | Vincent van Gogh, *The Starry Night* (1889). Public domain. [Wikimedia Commons](https://commons.wikimedia.org/wiki/File:Van_Gogh_-_Starry_Night_-_Google_Art_Project.jpg), 960 px rendition resized to 800x633 (bilinear) and re-encoded as JPEG quality 85 with stb_image_write. |

The visualiser orders images from easiest to hardest by measurement, not by this table: the average energy of the pixels that the DP carve removes (see `tools/make_viz.py`).

Do not add copyrighted artwork (for example, *The Persistence of Memory* is not in the public domain).
