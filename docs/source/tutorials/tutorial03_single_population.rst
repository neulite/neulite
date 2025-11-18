===================================================
Tutorial 3: 単一ポピュレーションネットワーク
===================================================

概要
====

このチュートリアルでは、複数のニューロンからなる単一ポピュレーションのシミュレーションを行います。

.. warning::
   bionet_liteは現在、外部かファイルを用いたスパイク入力には対応していません。このチュートリアルでは、生物物理学的ニューロンのポピュレーションファイル・シナプス接続情報ファイルが生成されることを確認します。


前提条件: BMTKチュートリアルの実行
==================================

.. important::
   **このチュートリアルを始める前に、必ずBMTKの対応するチュートリアルを先に実行してください**

   BMTKのチュートリアルを実行することで、必要なデータファイル（SWCファイル、JSONファイル等）が自動的にダウンロードされます。bionet_liteの実行には、BMTKのコード・データファイルが必要です。

**実行手順：**

1. **BMTKチュートリアルを実行**

   まず、以下のBMTKチュートリアルを最後まで実行してください：

   `Tutorial: Multi-Cell, Single Population Network (with BioNet) <https://alleninstitute.github.io/bmtk/tutorials/tutorial_03_single_pop.html>`_


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

生成されるファイル
==================

bionet_liteは以下のディレクトリ内にファイルを生成します：

* スクリプトに ``base_dir`` が設定されている場合: ``{base_dir}_nl/``
* ``base_dir`` が設定されていない場合: ``neulite/``

**ディレクトリ構造**

.. code-block:: text

   neulite/  (または {base_dir}_nl/)
   ├── mcortex_population.csv
   ├── mcortex_mcortex_connection.csv
   ├── kernel/
   │   └── config.h
   └── data/
       ├── Scnn1a_473845048_m.swc
       └── 472363762_fit.csv

mcortex_population.csv
----------------------

100個のニューロンの情報が含まれます。カラム構成は以下の通りです：

.. code-block:: text

   #n_cell,n_comp,name,swc_file,ion_file
   100,3682,Scnn1a_100,data/Scnn1a_473845048_m.swc,data/472363762_fit.csv

* ``#n_cell``: 細胞数（100個）
* ``n_comp``: コンパートメント数
* ``name``: ポピュレーション名
* ``swc_file``: SWCモーフォロジファイルのパス
* ``ion_file``: イオンチャネルパラメータファイルのパス

mcortex_mcortex_connection.csv
-------------------------------

シナプス接続情報が含まれます。bionet_liteは、bionetが定義した接続ルールを適用し、具体的な接続位置（post cid）を決定します。

.. code-block:: text

   #pre nid,post nid,post cid,weight,tau_decay,tau_rise,erev,delay,e/i
   0,5,1523,0.00005,1.7,0.1,0.0,2,e
   0,12,892,0.00005,1.7,0.1,0.0,2,e
   ...

* ``#pre nid``: シナプス前ニューロンID
* ``post nid``: シナプス後ニューロンID
* ``post cid``: シナプス後コンパートメントID
* ``weight``: シナプス重み
* ``tau_decay``: 減衰時定数（ms）
* ``tau_rise``: 上昇時定数（ms）
* ``erev``: 逆転電位（mV）
* ``delay``: シナプス遅延（ms）
* ``e/i``: 興奮性(e)または抑制性(i)

仕様
----------

bionet_liteは、save_edges()が呼ばれたタイミングでシナプス接続位置を決定します。これは通常bionet.simulatorが実行時に行う処理です。事前に定められた接続ルールに従うコンパートメントからランダムに接続を作成します。


次のステップ
============

* :doc:`tutorial04_multiple_populations` - 複数ポピュレーションのネットワーク
* :doc:`../user_guide/basic_usage` - 基本的な使い方
