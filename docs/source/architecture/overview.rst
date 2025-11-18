==================
概要
==================

システム構成
============

bionet_liteとNeuliteは、BMTKのモジュラー設計を活用したシステムです。

全体像
------

.. image:: ../_static/763493463-bmtk.png
   :width: 700px
   :align: center

.. raw:: html

   <div style="margin-bottom: 30px;"></div>

BMTKは以下の3つの主要コンポーネントに分かれています：

* **Builder**: ネットワーク構築
* **Simulator**: シミュレーション実行
* **Analyzer**: 結果の解析

bionet_liteはBuilderモジュールを拡張し、Neuliteは新しいSimulatorとして機能します。

詳細度に応じたインターフェース
==============================

BMTKは、異なる詳細度のモデルに対応するため、複数のインターフェースを提供しています：

* **bionet**: 生物物理学的ネットワーク（`NEURON <https://www.neuron.yale.edu/>`_）
* **bionet_lite**: 生物物理学的ネットワーク（Neulite） **← 本プロジェクト**
* **pointnet**: ポイントニューロンネットワーク（NEST）
* **filternet**: フィルターネットワーク（LGN Model）
* **popnet**: ポピュレーション統計ネットワーク（DiPDE）

bionet_liteの位置づけ
======================

bionet_liteは、bionetと同じく生物物理学的ネットワークを扱いますが、以下の点で異なります：

* **ターゲットシミュレータ**: `NEURON <https://www.neuron.yale.edu/>`_ → Neulite
* **出力形式**: `SONATA <https://github.com/AllenInstitute/sonata/tree/master>`_ 形式 + Neulite専用CSV形式
* **接続位置の決定**: シミュレータ実行時 → ネットワーク構築時

設計哲学
========

ネットワーク構築とシミュレーションの分離
----------------------------------------

BMTKの設計哲学に従い、bionet_liteは以下の利点を持ちます：

1. **エラーハンドリングの向上**: ネットワーク構築時のエラーをシミュレーション前に検出
2. **再現性**: 同じネットワークファイルで何度でも同じシミュレーションを実行可能
3. **柔軟性**: ネットワークファイルを手動で編集・検証可能

軽量化の実現
------------

Neuliteカーネルの軽量化という設計指針に従い、bionet_liteは以下の処理を行います：

* **前処理の集約**: 形態ファイル、イオンチャネル設定の前処理をビルド時に実行
* **シンプルなファイル形式**: CSV形式による読み書きのしやすさ
* **明示的なデータ構造**: シミュレータ側での複雑な処理を不要に

次のステップ
============

* :doc:`design_and_implementation` - BMTKとBioNetの背景、bionet_liteの設計と実装詳細
