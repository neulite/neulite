===============================================
Tutorial 2: スパイク入力シミュレーション
===============================================

概要
====

このチュートリアルでは、生物物理学的ニューロンに対してスパイク入力を行う方法を説明します。

.. warning::
   bionet_liteは現在、Virtual cellによる入力には対応していません。このチュートリアルでは、生物物理学的ニューロンのポピュレーション情報ファイルが生成されることを確認します。

前提条件: BMTKチュートリアルの実行
==================================

.. important::
   **このチュートリアルを始める前に、必ずBMTKの対応するチュートリアルを先に実行してください**

   BMTKのチュートリアルを実行することで、必要なデータファイル（SWCファイル、JSONファイル等）が自動的にダウンロードされます。bionet_liteの実行には、BMTKのコード・データファイルが必要です。

**実行手順：**

1. **BMTKチュートリアルを実行**

   まず、以下のBMTKチュートリアルを最後まで実行してください：

   `Tutorial: Single Cell Simulation with External Feedforward Input (with BioNet) <https://alleninstitute.github.io/bmtk/tutorials/tutorial_02_single_cell_syn.html>`_

2. **bionet_liteのコード例を実行**

   BMTKチュートリアルが完了したら、次のセクションに進んでください。

BMTKからbionet_liteへの変更
==============================

**変更が必要なのはimport文だけです。** それ以外のコードは完全に同じです。

**変更前（BMTK）:**

.. code-block:: python

   from bmtk.builder.networks import NetworkBuilder

**変更後（bionet_lite）:**

.. code-block:: python

   from bionetlite import NeuliteBuilder as NetworkBuilder

コード例
==============

.. code-block:: python

   from bionetlite import NeuliteBuilder as NetworkBuilder

   # Create network
   net = NetworkBuilder('mcortex')

   # Add biophysical neurons
   net.add_nodes(
       N=1,
       pop_name='Scnn1a',
       model_type='biophysical',
       model_template='ctdb:Biophys1.hoc',
       dynamics_params='472363762_fit.json',
       morphology='Scnn1a_473845048_m.swc'
   )

   # Add Virtual cells for input
   # Note: bionet_lite ignores Virtual cells
   net.add_nodes(
       N=10,
       pop_name='input',
       model_type='virtual'
   )

   # Build network
   net.build()

   # Save to files
   net.save_nodes(output_dir='sim_ch02/network')

生成されるファイル
==================

このシミュレーションでは、生物物理学的ニューロンのpopulationファイルのみが生成されます。Virtual cellの情報は含まれません。

mcortex_population.csv
----------------------

生物物理学的ニューロン（1個）の情報が含まれます。カラム構成は以下の通りです：

.. code-block:: text

   #n_cell,n_comp,name,swc_file,ion_file
   1,3682,Scnn1a_100,data/Scnn1a_473845048_m.swc,data/472363762_fit.csv

* ``#n_cell``: Number of cells
* ``n_comp``: Number of compartments
* ``name``: Population name
* ``swc_file``: Path to SWC morphology file
* ``ion_file``: Path to ion channel parameter file

.. note::
   Virtual cellは出力されないため、``mcortex_mcortex_connection.csv`` はヘッダーのみのファイルとなります。

制限事項
========

bionet_liteならびにNeuliteは、以下の入力方式に対応しています：

* ✅ 定常電流入力
* ❌ ファイルからのスパイク列入力（現在非対応）
* ❌ Virtual cellによる入力（現在非対応）

.. note::
   スパイク入力機能が必要な場合は、Neulite側での実装が必要です。

次のステップ
============

* :doc:`tutorial03_single_population` - 複数ニューロンを含む単一ポピュレーション
* :doc:`../advanced/specification` - bionet_liteの詳細な仕様と制限事項
