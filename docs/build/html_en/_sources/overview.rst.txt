========
概要
========

このドキュメントについて
========================

このドキュメントでは、生物物理学的ニューロンモデルのための高性能シミュレーション環境である **Neulite** と、そのフロントエンドである **bionet_lite** の使用方法を説明します。


.. image:: _static/767369955-bionet_lite_overview.png
   :width: 700px
   :align: center

.. raw:: html

   <div style="margin-bottom: 30px;"></div>



Neuliteは、`Brain Modeling Toolkit (BMTK) <https://github.com/AllenInstitute/bmtk>`_ で記述されたモデルを **少量の書き換えで動作可能** にします。既存のBMTKモデルを活用しながら、高性能な計算環境での実行を実現します。



システム構成
============

Neuliteシステムは、フロントエンドとバックエンドが明確に分離された設計になっています。この分離により、Pythonでの柔軟なネットワーク構築とC言語による高速なシミュレーション実行を両立しています。

* **バックエンド**: Neuliteカーネル（Cベース）
* **フロントエンド**: Bionet_liteモジュール（Pythonベース）

この分離により、ネットワーク構築と実際のシミュレーション実行が独立しており、柔軟な開発とメンテナンスが可能です。



.. image:: _static/763493463-bmtk.png
   :width: 700px
   :align: center

.. raw:: html

   <div style="margin-bottom: 30px;"></div>

Bionet_liteとは
===============

**Bionet_lite** は、`Brain Modeling Toolkit (BMTK) <https://github.com/AllenInstitute/bmtk>`_ と互換性を持つネットワーク構築ツールです。既存のBMTKコードを最小限の変更で再利用することで、`SONATA <https://github.com/AllenInstitute/sonata/tree/master>`_ などの標準データ形式を自動的にサポートします。

   * Pythonで実装されたネットワーク構築・前処理ツール
   * SWCファイルの前処理（軸索の修正、DFSソート）を実行
   * イオンチャネル設定をJSON形式からCSV形式に変換
   * Neuliteカーネル用の設定ファイルを自動生成

詳細は :doc:`bionet_lite` をご覧ください。

Neuliteカーネルとは
====================

**Neuliteカーネル** は、大規模な生物物理学的ニューロンモデルのシミュレーションを高速に実行するための軽量シミュレータです。`Allen Cell Types Database <https://celltypes.brain-map.org/>`_ の `Perisomatic model <https://brain-map.org/our-research/computational-modeling/perisomatic-biophysical-single-neuron-models>`_ に特化し、ラズベリーパイから「`富岳 <https://www.r-ccs.riken.jp/en/fugaku/>`_」まで幅広い環境で動作します。


   * C言語（C17）で記述された高速シミュレータ
   * bionet_liteが生成した設定ファイルを読み込み
   * ユーザーがC言語で機能を拡張可能


詳細は :doc:`neulite` をご覧ください。


使用の流れ
==========

**準備**: BMTKのネットワーク構築コードを用意する

1. **ネットワーク構築**: bionet_liteをimportしたPythonコードを実行（Neuliteカーネル用の設定ファイルが自動生成される）
2. **シミュレーション実行**: Neuliteカーネルでシミュレーションを実行
3. **結果の解析**: 生成されたデータを解析

次のステップ
============

* :doc:`neulite` - Neuliteシミュレータの詳細
* :doc:`bionet_lite` - Bionet_liteの詳細と背景
* :doc:`getting_started/index` - セットアップとクイックスタート
* :doc:`tutorials/tutorial01_single_cell` - 実践的なチュートリアル
* :doc:`architecture/overview` - 詳細なアーキテクチャの説明
