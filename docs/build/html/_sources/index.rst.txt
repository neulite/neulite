================================================
Neulite プロジェクトドキュメンテーション
================================================

プロジェクトの概要
==================

Neuliteプロジェクトは、スーパーコンピュータ「`富岳 <https://www.r-ccs.riken.jp/en/fugaku/>`_」で大規模な生物物理学的ニューロンシミュレーションを実行するために開発されたシステムです。**900万個のニューロンと2600億個のシナプス** を含むマウス全皮質モデルのシミュレーションを実現しています（詳細は :doc:`publication` を参照してください）。


.. image:: _static/teaser.png
   :width: 700px
   :align: center

.. raw:: html

   <div style="margin-bottom: 30px;"></div>

このドキュメントでは、生物物理学的ニューロンモデルのための高性能シミュレーション環境である **Neulite** と、そのフロントエンドである **bionet_lite** について説明します。


BMTKとNeuliteの関係
====================

Neuliteは、`Brain Modeling Toolkit (BMTK) <https://github.com/AllenInstitute/bmtk>`_ で記述されたモデルを **少量の書き換えで動作可能** にします。既存のBMTKモデルを活用しながら、高性能な計算環境での実行を実現します。


主な特徴
========

* 🚀 **既存コードの再利用**: BMTKで書かれたコードのimport文を変更するだけで利用可能
* 📊 **データ標準化**: `SONATA <https://github.com/AllenInstitute/sonata/tree/master>`_ や `Allen Cell Types Database <https://celltypes.brain-map.org/>`_ などの標準形式を自動サポート
* 🔧 **軽量設計**: UNIX哲学に基づくシンプルで理解しやすいカーネル実装
* 🌐 **高い可搬性**: ラズベリーパイからスーパーコンピュータまで様々な環境で実行可能


コンポーネント
--------------



.. image:: _static/763493463-bmtk.png
   :width: 700px
   :align: center

.. raw:: html

   <div style="margin-bottom: 30px;"></div>


**Bionet_lite（フロントエンド）**
   `Brain Modeling Toolkit (BMTK) <https://github.com/AllenInstitute/bmtk>`_ のbionetモジュールを拡張したPythonパッケージ。既存のBMTKコードを少量の変更で利用可能。

**Neuliteカーネル（シミュレータ）**
   軽量で高性能な生物物理学的ニューロンシミュレータ。C言語で実装されたクリーンなコードベースであり、軽量で拡張性の高い設計が特徴。


クイックスタート
================

.. code-block:: python

   # bionetの代わりにbionet_liteをインポート
   from bionetlite import NeuliteBuilder as NetworkBuilder

   # あとは通常のbionetと同じように使用
   net = NetworkBuilder('v1')
   net.add_nodes(...)
   net.build()
   net.save_nodes(output_dir='network')
   net.save_edges(output_dir='network')

.. note::
   完全なセットアップ手順については、:doc:`getting_started/index` を参照してください。

ドキュメント目次
================

.. toctree::
   :maxdepth: 2
   :caption: はじめに

   overview
   bionet_lite
   neulite
   getting_started/index

.. toctree::
   :maxdepth: 2
   :caption: チュートリアル

   tutorials/tutorial01_single_cell
   tutorials/tutorial02_spike_input
   tutorials/tutorial03_single_population
   tutorials/tutorial04_multiple_populations

.. toctree::
   :maxdepth: 2
   :caption: ユーザーガイド

   user_guide/basic_usage
   user_guide/configuration
   user_guide/parallel_execution

.. toctree::
   :maxdepth: 2
   :caption: アーキテクチャ

   architecture/overview
   architecture/design_and_implementation

.. toctree::
   :maxdepth: 2
   :caption: APIリファレンス

   api_reference/index
   api_reference/file_formats

.. toctree::
   :maxdepth: 2
   :caption: 高度なトピック

   advanced/allen_v1_model
   advanced/specification

.. toctree::
   :maxdepth: 1
   :caption: その他

   faq
   publication
   how_to_cite
   changelog

インデックスと検索
==================

* :ref:`genindex`
* :ref:`search`
