==============
ファイル形式
==============

bionet_liteが生成するファイルの形式について説明します。

概要
====

bionet_liteは、以下の2種類のファイルを生成します：

1. `SONATA <https://github.com/AllenInstitute/sonata/tree/master>`_ 形式: BMTKと互換性のあるHDF5形式
2. Neulite形式: Neuliteカーネルが読み込むCSV形式

SONATA形式ファイル
==================

BMTKの他のツールとの互換性のため、SONATA形式のファイルも生成されます。

**ノードファイル:**

* ``<network_name>_nodes.h5`` - HDF5形式のノードデータ
* ``<network_name>_node_types.csv`` - CSVノードタイプ定義

**エッジファイル:**

* ``<src>_<trg>_edges.h5`` - HDF5形式のエッジデータ
* ``<src>_<trg>_edge_types.csv`` - CSVエッジタイプ定義

これらのファイルはBMTKの標準形式に従います。

Neulite形式ファイル
====================

Neuliteカーネル用のファイルです。

<network_name>_population.csv
==============================

ニューロンポピュレーション情報を定義します。

ヘッダー
--------

.. code-block:: text

   #n_cell,n_comp,name,swc_file,ion_file

カラム
------

.. list-table::
   :header-rows: 1
   :widths: 20 15 65

   * - カラム名
     - 型
     - 説明
   * - #n_cell
     - int
     - このモデルタイプのニューロン数
   * - n_comp
     - int
     - コンパートメント数（SWCファイルの行数）
   * - name
     - string
     - モデル名（形式: ``<pop_name>_<node_type_id>``）
   * - swc_file
     - string
     - 形態ファイルのパス（例: ``data/Scnn1a.swc``）
   * - ion_file
     - string
     - イオンチャネル設定ファイルのパス（例: ``data/472363762_fit.csv``）

例
--

.. code-block:: text

   #n_cell,n_comp,name,swc_file,ion_file
   80,195,Scnn1a_100,data/Scnn1a_473845048_m.swc,data/472363762_fit.csv
   20,158,PV_101,data/Pvalb_473862421_m.swc,data/472912177_fit.csv

<src>_<trg>_connection.csv
===========================

シナプス接続情報を定義します。

ヘッダー
--------

.. code-block:: text

   #pre nid,post nid,post cid,weight,tau_decay,tau_rise,erev,delay,e/i

カラム
------

.. list-table::
   :header-rows: 1
   :widths: 20 15 65

   * - カラム名
     - 型
     - 説明
   * - #pre nid
     - int
     - 送信元ニューロンID（0から始まる）
   * - post nid
     - int
     - 受信先ニューロンID（0から始まる）
   * - post cid
     - int
     - 受信先コンパートメントID（SWCファイルのノード番号）
   * - weight
     - float
     - シナプス重み
   * - tau_decay
     - float
     - 減衰時定数（ms）、exp2synの ``tau2``
   * - tau_rise
     - float
     - 立ち上がり時定数（ms）、exp2synの ``tau1``
   * - erev
     - float
     - 平衡電位（mV）
   * - delay
     - int
     - シナプス遅延（ms、整数に丸められる）
   * - e/i
     - string
     - 'e'（興奮性）または'i'（抑制性）

.. note::
   * ``post cid`` は ``target_sections`` で指定されたセクションからランダムにサンプリングされます
   * ``delay`` は整数値にROUND_HALF_UPで丸められます（2.5は3になる）

例
--

.. code-block:: text

   #pre nid,post nid,post cid,weight,tau_decay,tau_rise,erev,delay,e/i
   0,5,42,0.0005,1.7,0.1,0.0,2,e
   0,8,15,0.0005,1.7,0.1,0.0,2,e
   25,5,12,0.002,8.3,0.5,-70.0,1,i

swcファイル（処理済み）
=======================

Neulite用に前処理されたニューロン形態ファイルです。

処理内容
--------

bionet_liteは以下の処理を行います：

1. **ゼロベースインデックスへの変換** - 1ベースを0ベースに変換
2. **軸索の修正** - `Perisomatic model <https://brain-map.org/our-research/computational-modeling/perisomatic-biophysical-single-neuron-models>`_ への変換（簡略化された人工軸索を追加）
3. **深さ優先探索（DFS）ソート** - Neuliteでの効率的な処理のため

形式
----

.. code-block:: text

   #id type x y z r parent
   0 1 0.0 0.0 0.0 12.5 -1
   1 1 0.0 0.0 5.0 12.5 0
   2 2 0.0 0.0 35.0 0.5 1
   3 2 0.0 0.0 65.0 0.5 2
   4 3 -5.0 0.0 5.0 0.8 1

**カラム:**

* ``#id`` - コンパートメント番号（0から始まる）
* ``type`` - 1=soma, 2=axon, 3=basal dendrite, 4=apical dendrite
* ``x, y, z`` - 3D座標（μm）
* ``r`` - 半径（μm）
* ``parent`` - 親のコンパートメント番号（-1はルート）

出力先: ``<neulite_dir>/data/``

イオンチャネル設定CSV
======================

イオンチャネルパラメータファイル（JSON→CSV変換済み）です。

ヘッダー
--------

.. code-block:: text

   Cm,Ra,leak,e_pas,gamma,decay,NaV,NaTs,NaTa,Nap,Kv2like,Kv3_1,K_P,K_T,Kd,Im,Im_v2,Ih,SK,Ca_HVA,Ca_LVA

行の意味
--------

4行のCSVファイルで、各行はニューロンのセクションを表します：

.. code-block:: text

   1,<soma values>
   2,<axon values>
   3,<basal dendrite values>
   4,<apical dendrite values>

主要パラメータ
--------------

* **Cm** (μF/cm²) - 膜容量
* **Ra** (Ω·cm) - 軸索抵抗
* **leak** (S/cm²) - リークコンダクタンス（g_pas）
* **e_pas** (mV) - リーク平衡電位
* **gamma, decay** - カルシウムダイナミクス
* **イオンチャネル** - NaV, NaTs, NaTa, Nap, Kv2like, Kv3_1, K_P, K_T, Kd, Im, Im_v2, Ih, SK, Ca_HVA, Ca_LVA

出力先: ``<neulite_dir>/data/``

config.h
========

Neuliteカーネル用の設定ヘッダーファイルです。

生成されるマクロ
----------------

.. code-block:: c

   #pragma once

   // Simulation parameters
   #define TSTOP ( 3000.0 )
   #define DT ( 0.1 )
   #define INV_DT ( ( int ) ( 1.0 / ( DT ) ) )

   // Neuron parameters
   #define SPIKE_THRESHOLD ( -15.0 )
   #define ALLACTIVE ( 0 )

   // Current injection parameters
   #define I_AMP ( 0.1 )
   #define I_DELAY ( 500.0 )
   #define I_DURATION ( 500.0 )

パラメータ
----------

* **TSTOP** (ms) - シミュレーション終了時間
* **DT** (ms) - タイムステップ
* **INV_DT** - タイムステップの逆数（整数）
* **SPIKE_THRESHOLD** (mV) - スパイク検出閾値
* **ALLACTIVE** - allactiveモデルフラグ（0 or 1）
* **I_AMP** (nA) - 電流注入の振幅
* **I_DELAY** (ms) - 電流注入の開始時刻
* **I_DURATION** (ms) - 電流注入の持続時間

出力先: ``<neulite_dir>/kernel/config.h``

ディレクトリ構造
================

bionet_liteは以下のディレクトリ構造を生成します：

.. code-block:: text

   project/
   ├── network/                      # SONATA形式（BMTK互換）
   │   ├── V1_nodes.h5
   │   ├── V1_node_types.csv
   │   ├── V1_V1_edges.h5
   │   └── V1_V1_edge_types.csv
   └── neulite/                      # Neulite形式
       ├── V1_population.csv
       ├── V1_V1_connection.csv
       ├── kernel/
       │   └── config.h
       └── data/
           ├── *.swc                 # 変換されたSWCファイル
           └── *.csv                 # イオンチャネルCSV

次のステップ
============

* :doc:`index` - APIリファレンス
* :doc:`../user_guide/basic_usage` - 基本的な使用方法
* :doc:`../user_guide/configuration` - 設定ファイルの詳細
