============
Publication
============

このページでは、Neuliteとbionet_liteに関する主要な論文を紹介します。

主要論文
========
1
Kuriyama, R.*, Akira, K.*, Green, L., Herrera, B., Dai, K., Iura, M., Gouaillardet, G., Terasawa, A., Kobayashi, T., Igarashi, Jun., Arkhipov, A., & Yamazaki, T. (2025).
**Microscopic-Level Mouse Whole Cortex Simulation Composed of 9 Million Biophysical Neurons and 26 Billion Synapses on the Supercomputer Fugaku.**
In *Proceedings of the International Conference for High Performance Computing, Networking, Storage and Analysis (SC '25)*, November 16–21, 2025, St Louis, MO, USA. ACM, New York, NY, USA, 11 pages.

DOI: `10.1145/3712285.3759819 <https://doi.org/10.1145/3712285.3759819>`_

\* 共同第一著者

論文の概要
==========

この論文では、スーパーコンピュータ「`富岳 <https://www.r-ccs.riken.jp/en/fugaku/>`_」を使用して、マウス全皮質の顕微鏡レベルでのシミュレーションを実現した研究を報告しています。

主な成果
--------

* **900万個の生物物理学的ニューロン** を含むマウス全皮質モデル
* **2600億個のシナプス** による複雑な接続構造
* 「富岳」での効率的な並列実行
* `Perisomatic model <https://brain-map.org/our-research/computational-modeling/perisomatic-biophysical-single-neuron-models>`_ を用いた詳細な生物物理学的シミュレーション

技術的貢献
----------

* **bionet_lite**: BMTKを拡張したネットワーク構築フレームワーク
* **Neulite**: 軽量で高性能な生物物理学的ニューロンシミュレータ
* **最適化技術**: SVE（Scalable Vector Extensions）を活用したSIMD計算
