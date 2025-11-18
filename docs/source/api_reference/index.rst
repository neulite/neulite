================
APIリファレンス
================

このセクションでは、bionet_liteの主要なAPIについて説明します。

概要
====

bionet_liteは、BMTKのネットワークビルダーを拡張したPythonライブラリです。
BMTKコードのimport文を変更するだけで、Neuliteシミュレータ用のファイルを生成できます。

NeuliteBuilderクラス
====================

bionet_liteのメインインターフェースです。BMTKの ``NetworkBuilder`` を継承しています。

初期化
------

.. code-block:: python

   from bionetlite import NeuliteBuilder

   net = NeuliteBuilder(
       name,                        # Network name (required)
       convert_morphologies=True,   # Convert SWC files for Neulite
       convert_ion_channels=True,   # Convert ion channel JSON to CSV
       simulation_config=None,      # Simulation config (optional)
       generate_config_h=True       # Auto-generate config.h
   )

主要メソッド
------------

add_nodes()
^^^^^^^^^^^

ネットワークにノード（ニューロン）を追加します。

.. code-block:: python

   net.add_nodes(
       N=100,                              # Number of neurons
       pop_name='Exc',                     # Population name
       model_type='biophysical',           # Model type
       morphology='Scnn1a.swc',            # Morphology file name
       dynamics_params='472363762_fit.json',  # Ion channel parameters
       ei='e'                              # 'e' (excitatory) or 'i' (inhibitory)
   )

**対応モデル**: 現在は ``model_type='biophysical'`` のみ対応

add_edges()
^^^^^^^^^^^

ネットワークにエッジ（シナプス接続）を追加します。

.. code-block:: python

   net.add_edges(
       source={'pop_name': 'Exc'},         # Source population
       target={'pop_name': 'Inh'},         # Target population
       connection_rule=5,                  # Number of connections or probability
       syn_weight=0.001,                   # Synaptic weight
       delay=2.0,                          # Delay (ms, rounded to integer)
       target_sections=['somatic', 'basal'],  # Target sections
       dynamics_params='ExcToInh.json'     # Synapse parameters
   )

**対応モデル**: 現在は ``exp2syn`` のみ対応

build()
^^^^^^^

ネットワークをビルドします。

.. code-block:: python

   net.build()

save_nodes()
^^^^^^^^^^^^

ノード情報を保存します（`SONATA <https://github.com/AllenInstitute/sonata/tree/master>`_ 形式 + Neulite形式）。

.. code-block:: python

   net.save_nodes(output_dir='network')

**生成されるファイル**:

* SONATA形式: ``<name>_nodes.h5``, ``<name>_node_types.csv``
* Neulite形式: ``<name>_population.csv``, SWCファイル、イオンチャネルCSV

save_edges()
^^^^^^^^^^^^

エッジ情報を保存します（SONATA形式 + Neulite形式）。

.. code-block:: python

   net.save_edges(output_dir='network')

**生成されるファイル**:

* SONATA形式: ``<src>_<trg>_edges.h5``, ``<src>_<trg>_edge_types.csv``
* Neulite形式: ``<src>_<trg>_connection.csv``

設定メソッド
------------

set_config_path()
^^^^^^^^^^^^^^^^^

config.jsonのパスを設定します。

.. code-block:: python

   net.set_config_path('config.json')

set_dir()
^^^^^^^^^

Neuliteファイルのディレクトリパスを設定します。

.. code-block:: python

   net.set_dir(ion_dir='data/ion', swc_dir='data/morphology')

基本的な使用例
==============

.. code-block:: python

   from bionetlite import NeuliteBuilder

   # Network initialization
   net = NeuliteBuilder('V1')

   # Add excitatory neurons
   net.add_nodes(N=100, pop_name='Exc', model_type='biophysical',
                 morphology='Scnn1a.swc',
                 dynamics_params='472363762_fit.json', ei='e')

   # Add inhibitory neurons
   net.add_nodes(N=25, pop_name='Inh', model_type='biophysical',
                 morphology='Pvalb.swc',
                 dynamics_params='472912177_fit.json', ei='i')

   # Exc -> Exc connection
   net.add_edges(source={'pop_name': 'Exc'}, target={'pop_name': 'Exc'},
                 connection_rule=5, syn_weight=0.0005, delay=2.0,
                 target_sections=['basal', 'apical'],
                 dynamics_params='ExcToExc.json')

   # Exc -> Inh connection
   net.add_edges(source={'pop_name': 'Exc'}, target={'pop_name': 'Inh'},
                 connection_rule=3, syn_weight=0.001, delay=2.0,
                 target_sections=['somatic', 'basal'],
                 dynamics_params='ExcToInh.json')

   # Build and save
   net.build()
   net.save_nodes(output_dir='network')
   net.save_edges(output_dir='network')

並列実行
========

大規模ネットワークの構築は、MPIで並列実行できます。

.. code-block:: bash

   mpirun -n 4 python build_network.py

次のステップ
============

* :doc:`file_formats` - 生成されるファイル形式の詳細
* :doc:`../user_guide/basic_usage` - 実践的な使用例
* :doc:`../tutorials/tutorial01_single_cell` - チュートリアル
* :doc:`../advanced/specification` - 仕様と制限事項
