===================
Neuliteカーネルとは
===================

概要
====

Neuliteカーネルは、生物物理学的ニューロンモデルとネットワークモデルのための **軽量なシミュレータ** です。`Allen Cell Types Database <https://celltypes.brain-map.org/>`_ のニューロンモデルのシミュレーションに特化しています。

特徴
====


軽量設計
----------

NeuliteはC17で記述されたクリーンなコードベースを特徴としており、軽量かつ高い拡張性を持ちます。


高い可搬性
----------

Neuliteは ラズベリーパイからスーパーコンピュータまで、様々なアーキテクチャで動作します。この高い可搬性により、ネットワーク構築とシミュレーション実行環境を分離することが可能です。
例えば、依存ライブラリの多いネットワーク構築部分（Bionet_lite）はローカルで実行し、生成されたファイルを元にNeuliteカーネルのみをスーパーコンピュータ上で利用することができます。

技術的詳細
==========


設計思想
--------

Neuliteの設計は **UNIX哲学** の影響を受けており、以下の原則に従っています：

* **単純さの重視**: 複雑さを避け、理解しやすいコードを目指す
* **モジュール化**: 各コンポーネントが独立した機能を持つ
* **移植性**: 特定の環境に依存しない設計

対象ニューロンモデル
--------------------

Neuliteは `Perisomatic model <https://brain-map.org/our-research/computational-modeling/perisomatic-biophysical-single-neuron-models>`_ を採用したAllen Cell Types Databaseのニューロンモデルに特化しています。このモデルの詳細については、:doc:`advanced/specification` を参照してください。

使用環境
========

推奨環境
--------

* **環境**: Linux、macOSなど
* **コンパイラ**: C17対応のCコンパイラ（gcc 7.0以上推奨）


富岳での実績
------------

Neuliteは「`富岳 <https://www.r-ccs.riken.jp/en/fugaku/>`_」スーパーコンピュータ向けに最適化され、900万個の生物物理学的ニューロンと260億個のシナプスを含むマウス全皮質の顕微鏡レベルでのシミュレーションを実現した実績があります（:doc:`publication` を参照）。

**採用技術**:

* **Scalable Vector Extensions (SVE)**: SIMDベクトル計算を全面的に活用
* **並列実装**: 大規模ネットワークの効率的なシミュレーション
* **MPI対応**: ノード間並列処理

関連リンク
==========

* 公式サイト: https://numericalbrain.org/neulite/
* 引用情報: :doc:`how_to_cite`


.. seealso::
   Neuliteの技術的な詳細や制約については、:doc:`advanced/specification` を参照してください。
