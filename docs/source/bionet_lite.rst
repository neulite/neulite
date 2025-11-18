==================
Bionet_liteとは
==================

概要
====

Bionet_liteは、`Brain Modeling Toolkit (BMTK) <https://github.com/AllenInstitute/bmtk>`_ のbionetモジュールを拡張したPythonパッケージです。Neuliteカーネル用のネットワークデータを生成し、既存のBMTKコードを最小限の変更で再利用できます。

特徴
====

既存コードの再利用
------------------

BMTKで記述されたネットワーク構築コードのimport文を変更するだけで、Neulite用のデータを生成できます。

.. code-block:: python

   # BMTK case
   from bmtk.builder.networks import NetworkBuilder

   net = NetworkBuilder('v1')
   net.add_nodes(...)
   net.build()
   net.save_nodes(output_dir='network')

.. code-block:: python

   # bionet_lite case (only change import)
   from bionetlite import NeuliteBuilder as NetworkBuilder

   net = NetworkBuilder('v1')
   net.add_nodes(...)  # Same code
   net.build()         # Same code
   net.save_nodes(output_dir='network')  # Same code


技術的詳細
==========

BMTKのNetworkBuilderクラスを拡張し、bionetの全機能を維持しながらNeulite用の出力を追加します。

.. code-block:: python

   # bionet case
   from bmtk.builder.networks import NetworkBuilder

   # bionet_lite case
   from bionetlite import NeuliteBuilder as NetworkBuilder



主な処理
--------

* **SWCファイルの前処理**: `Perisomatic model <https://brain-map.org/our-research/computational-modeling/perisomatic-biophysical-single-neuron-models>`_ への変換、DFSソート
* **イオンチャネル設定の変換**: JSON形式からCSV形式へ
* **ネットワークファイルの生成**: <network_name>_population.csv、<src>_<trg>_connection.csv

詳細は :doc:`architecture/design_and_implementation` を参照してください。

出力ファイルの比較
==================

bionetとbionet_liteでは、生成されるファイルが異なります。

.. list-table:: 出力ファイルの比較
   :header-rows: 1
   :widths: 30 35 35

   * - ファイルの種類
     - bionet（BMTK）
     - bionet_lite
   * - ネットワークデータ
     - | ``<network>_nodes.h5``
       | ``<network>_node_types.csv``
       | ``<src>_<trg>_edges.h5``
       | ``<src>_<trg>_edge_types.csv``
     - | ``<network>_nodes.h5``
       | ``<network>_node_types.csv``
       | ``<src>_<trg>_edges.h5``
       | ``<src>_<trg>_edge_types.csv``
       | **+** ``<network>_population.csv``
       | **+** ``<src>_<trg>_connection.csv``
   * - 形態ファイル
     - 元のSWCファイル（変更なし）
     - **処理済みSWCファイル** （Perisomatic model変換、DFSソート済み）
   * - イオンチャネル
     - JSON形式（``*_fit.json``）
     - | JSON形式（``*_fit.json``）
       | **+** CSV形式（Neulite用に変換）

bionet_liteは、BMTKの標準出力（H5/CSV）を完全に維持しながら、Neulite用の追加ファイル（<network_name>_population.csv、<src>_<trg>_connection.csv）を生成します。これにより、同じネットワークデータで両方のシミュレータを実行できます。



使用例
======

.. code-block:: python

   from bionetlite import NeuliteBuilder as NetworkBuilder

   # Create network
   net = NetworkBuilder('v1')

   # Add nodes and edges (same as BMTK)
   net.add_nodes(...)
   net.add_edges(...)

   # Build and save
   net.build()
   net.save_nodes(output_dir='network')
   net.save_edges(output_dir='network')

データ標準化の自動サポート
--------------------------

このアプローチによって、BMTKが対応する `SONATA <https://github.com/AllenInstitute/sonata/tree/master>`_ や `Allen Cell Types Database <https://celltypes.brain-map.org/>`_ などの標準形式を自動的にサポートします。


関連リンク
==========

* Neuliteシミュレータ: :doc:`neulite`
* セットアップガイド: :doc:`getting_started/bionet_lite_setup`
* チュートリアル: :doc:`tutorials/tutorial01_single_cell`

.. seealso::
   Bionet_liteの設計詳細や実装、BMTKとBioNetの背景については、:doc:`architecture/design_and_implementation` を参照してください。
