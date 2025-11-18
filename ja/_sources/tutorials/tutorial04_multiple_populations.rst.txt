====================================================================
Tutorial 4: 複数のポピュレーションからなるネットワーク
====================================================================

概要
====

このチュートリアルでは、複数のニューロンタイプからなる複数ポピュレーションのシミュレーションを行います。


前提条件: BMTKチュートリアルの実行
==================================

.. important::
   **このチュートリアルを始める前に、必ずBMTKの対応するチュートリアルを先に実行してください**

   BMTKのチュートリアルを実行することで、必要なデータファイル（SWCファイル、JSONファイル等）が自動的にダウンロードされます。bionet_liteの実行には、BMTKのコード・データファイルが必要です。

**実行手順：**

1. **BMTKチュートリアルを実行**

   まず、以下のBMTKチュートリアルを最後まで実行してください：

   `Tutorial: Multi-Population Recurrent Network (with BioNet) <https://alleninstitute.github.io/bmtk/tutorials/tutorial_04_multi_pop.html>`_

2. **bionet_liteのコード例を実行**

   BMTKチュートリアルが完了したら、次のセクションに進んでください。

bionetからbionet_liteへの変更
==============================

**変更が必要なのはimport文だけです。** それ以外のコードは完全に同一です。

**変更前（BMTK）:**

.. code-block:: python

   from bmtk.builder.networks import NetworkBuilder

**変更後（bionet_lite）:**

.. code-block:: python

   from bionetlite import NeuliteBuilder as NetworkBuilder



生成されるファイル
==================

population.csv
--------------

100個のニューロン（興奮性80個 + 抑制性20個）の情報が含まれます。

connection.csv
--------------

4種類のシナプス接続（E→E、E→I、I→E、I→I）の情報が含まれます。

Neuliteが生成するファイル
=========================

bionet_liteは以下のディレクトリ内にファイルを生成します：

**ディレクトリ構造**

.. code-block:: text

   sim_ch04_nl/
   ├── V1_population.csv
   ├── V1_V1_connection.csv
   ├── kernel/
   │   └── config.h
   └── data/
       ├── Scnn1a_473845048_m.swc
       ├── Pvalb_470522102_m.swc
       ├── 472363762_fit.csv
       └── 472912177_fit.csv

V1_population.csv
-----------------

Neulite用のポピュレーションファイル。複数のニューロンタイプが定義されています：

.. code-block:: text

   #n_cell,n_comp,name,swc_file,ion_file
   80,3682,Scnn1a_100,data/Scnn1a_473845048_m.swc,data/472363762_fit.csv
   20,1900,PV_101,data/Pvalb_470522102_m.swc,data/472912177_fit.csv

各行は以下の情報を含みます：

* ``#n_cell``: 細胞数（80個の興奮性、20個の抑制性）
* ``n_comp``: コンパートメント数（Scnn1a: 3682、PV: 1900）
* ``name``: ポピュレーション名
* ``swc_file``: SWCモーフォロジファイルのパス
* ``ion_file``: イオンチャネルパラメータファイルのパス

V1_V1_connection.csv
--------------------

Neulite用の接続ファイル。ポピュレーション内およびポピュレーション間の全シナプス接続を定義します：

.. code-block:: text

   #pre nid,post nid,post cid,weight,tau_decay,tau_rise,erev,delay,e/i
   0,1,5,5e-05,1.0,3.0,0.0,2,e
   0,1,3351,5e-05,1.0,3.0,0.0,2,e
   0,2,28,5e-05,1.0,3.0,0.0,2,e
   ...

各カラムの意味：

* ``#pre nid``: シナプス前ニューロンID
* ``post nid``: シナプス後ニューロンID
* ``post cid``: シナプス後コンパートメントID
* ``weight``: シナプス重み
* ``tau_decay``: 減衰時定数
* ``tau_rise``: 上昇時定数
* ``erev``: 逆転電位
* ``delay``: シナプス遅延（ms）
* ``e/i``: 興奮性(e)または抑制性(i)

Neuliteカーネルの実行
=====================

前提条件
--------

Neuliteカーネルをビルドしておく必要があります。まだビルドしていない場合は、:doc:`../getting_started/neulite_build` を参照してビルドしてください。

ファイルの配置
--------------

bionet_liteが生成した ``neulite/`` ディレクトリの中身を、Neuliteカーネルのビルドディレクトリにコピーします：

.. code-block:: bash

   # Copy neulite directory contents to Neulite kernel directory
   cp -r neulite/* /path/to/neulite_kernel/

配置後のNeuliteカーネルディレクトリ構造は以下のようになります：

.. code-block:: text

   neulite_kernel/
   ├── nl                          # Neuliteカーネル実行ファイル（ビルド済み）
   ├── V1_population.csv           # bionet_liteが生成（コピー）
   ├── V1_V1_connection.csv        # bionet_liteが生成（コピー）
   ├── kernel/                     # bionet_liteが生成（コピー）
   │   └── config.h
   ├── data/                       # bionet_liteが生成（コピー）
   │   ├── Scnn1a_473845048_m.swc
   │   ├── Pvalb_470522102_m.swc
   │   ├── 472363762_fit.csv
   │   └── 472912177_fit.csv
   └── ... (その他のNeuliteカーネルのソースファイル)

.. note::
   Neuliteカーネルは ``kernel/config.h`` の設定を使用してビルドする必要があります。
   config.hをコピーした後、カーネルを再ビルドしてください。
   詳細は :doc:`../getting_started/neulite_build` を参照してください。

実行
----

Neuliteカーネルのディレクトリに移動して、シミュレーションを実行します：

.. code-block:: bash

   # Move to Neulite kernel directory
   cd /path/to/neulite_kernel

   # Run Neulite kernel
   ./nl V1_population.csv V1_V1_connection.csv

シミュレーション結果ファイル
----------------------------

Neuliteカーネル実行後、以下のテキストファイルが生成されます：

**s.dat** - スパイク時刻データ

スペース区切りのテキスト形式で、各行に1つのスパイクイベントを記録：

.. code-block:: text

   時刻(ms) ニューロンID
   12.5 0
   15.3 15
   18.9 42
   23.1 0
   ...

* 1列目: スパイクが発生した時刻（ms）
* 2列目: スパイクを発生させたニューロンのID（0-79: 興奮性、80-99: 抑制性）

**v.dat** - 膜電位データ

スペース区切りのテキスト形式で、各行に全ニューロンの膜電位を記録：

.. code-block:: text

   時刻(ms) ニューロン0 ニューロン1 ... ニューロン99
   0.0 -65.0 -65.0 ... -65.0
   0.1 -64.8 -64.9 ... -65.1
   0.2 -64.6 -64.8 ... -65.0
   ...

* 1列目: 時刻（ms）
* 2列目以降: 各ニューロンのソーマ膜電位（mV）、全100個のニューロン分

.. note::
   興奮性ニューロン（0-79）と抑制性ニューロン（80-99）のデータは同じファイル内に含まれています。

LIFモデルとの混在について
=========================

.. note::
   Tutorial 4には元々LIFモデルのニューロンも含まれますが、bionet_liteならびにNeuliteは現在biophysicalニューロンのみに対応しています。

   bionet_liteは以下のニューロンタイプに対応しています：

   * ✅ biophysical neurons
   * ❌ point neurons (LIFなど)
   * ❌ Virtual neurons

   そのため、このチュートリアルではbiophysicalニューロンに関連するファイルのみが生成されます。

次のステップ
============

* :doc:`../advanced/allen_v1_model` - 実際の研究モデルでの応用例
* :doc:`../user_guide/basic_usage` - 基本的な使い方
* :doc:`../advanced/specification` - サポートされる機能の詳細
