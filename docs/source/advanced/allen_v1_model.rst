==================================
Allen V1 Model: 大規模ネットワーク
==================================

概要
====

このチュートリアルでは、マウス一次視覚野カラムモデルをbionet_liteで構築する実例を示します。これは実際の研究で使用される大規模な生物物理学的モデルです。

目的
====

* 実際の研究モデルでのbionet_liteの使用方法を学ぶ
* 大規模ネットワークの構築プロセスを理解する
* Neuliteでのシミュレーション実行方法を学ぶ

.. figure:: ../_static/v1_stride_50_transparent.png
   :width: 600px
   :align: center

   Allen V1モデルのネットワーク構造。複数の皮質層にわたる生物物理学的ニューロンとそれらの接続を示す。

モデルの概要
============

`Allen V1モデル <https://brain-map.org/our-research/computational-modeling/models-of-the-mouse-primary-visual-cortex>`_ は、BMTKとbionetを用いて実装された生物物理学的モデルです。このモデルには以下の特徴があります：

* 形状を持った生物物理学的ニューロン
* 形状を持たないポイントニューロン（LIF等）
* 複数の皮質層を持つ階層構造
* 複雑なシナプス接続パターン

モデルの出典
------------

このモデルは、以下の論文で発表されたマウス一次視覚野の統合モデルに基づいています：

Billeh, Y. N., Cai, B., Gratiy, S. L., Dai, K., Iyer, R., Gouwens, N. W., Abbasi-Asl, R., Jia, X., Siegle, J. H., Olsen, S. R., Koch, C., Mihalas, S., & Arkhipov, A. (2020). **Systematic Integration of Structural and Functional Data into Multi-scale Models of Mouse Primary Visual Cortex.** *Neuron*, 106(3), 388-403.e18. https://doi.org/10.1016/j.neuron.2020.01.040

この論文では、構造データと機能データを体系的に統合した、マウス一次視覚野の生物物理学的モデルのネットワークが提示されています。

bionet_liteでの対応範囲
-----------------------

bionet_liteならびにNeuliteはbiophysicalモデルのみに対応しているため、このモデルに含まれる生物物理学的ニューロンからなるネットワークの構築に焦点を当てます。

.. note::
   元のモデルには約52000個の生物物理学的ニューロンが含まれています。

コード例
========

既存のBMTKコードの変更
----------------------

既存のAllen V1モデルのネットワーク構築コードでbionet_liteを使用するには、import文を変更するだけです：

.. code-block:: python

   # Before change
   # from bmtk.builder.networks import NetworkBuilder

   # After change
   from bionetlite import NeuliteBuilder as NetworkBuilder

   # Following code can be used as-is
   net = NetworkBuilder('v1')
   # ... network construction code ...

並列実行
========

大規模なネットワークの構築は計算量と時間がかかります。bionet_liteはbionetと同様に並列実行に対応しています：

.. code-block:: bash

   # Parallel execution with 4 processes
   mpirun -n 4 python build_network.py

生成されるファイル
==================

このモデルでは以下のファイルが生成されます：

* ``V1_population.csv`` - 約52000個のbiophysicalニューロンの情報
* ``V1_V1_connection.csv`` - 約1600万のシナプス接続情報
* 処理済みswcファイル群
* イオンチャネル設定ファイル群

結果の確認
==========

.. note::
   このセクションには、bionet_lite+Neuliteでの実行結果と、bionet+NEURONでの実行結果の比較図を追加してください。

   以下の図を含めることを推奨します：

   * ラスタープロット
   * 層ごとの発火率比較
   * 実行時間の比較

次のステップ
============

* :doc:`../user_guide/parallel_execution` - 並列実行の詳細



関連リンク
----------

* Allen Institute V1モデル: https://brain-map.org/our-research/computational-modeling/models-of-the-mouse-primary-visual-cortex
* Allen Cell Types Database: https://celltypes.brain-map.org/

.. warning::
  完全なコードとデータは、元のAllen Instituteのモデルリポジトリを参照してください。
