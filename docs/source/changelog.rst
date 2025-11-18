=========
変更履歴
=========

バージョン 1.0.0
================

*初回リリース*

Neuliteプロジェクトの初回公開版です。スーパーコンピュータ「`富岳 <https://www.r-ccs.riken.jp/en/fugaku/>`_」での大規模生物物理学的ニューロンシミュレーションを実現します。

主な成果
--------

**大規模シミュレーション実績**

* 900万個のニューロンと2600億個のシナプスを含むマウス全皮質モデルのシミュレーションを実現
* 富岳での実行に最適化

**bionet_lite**

ネットワーク構築ツール：

* BMTKの `NetworkBuilder` を拡張した `NeuliteBuilder` クラス
* 既存のBMTKコードをimport文の変更だけで利用可能
* biophysicalニューロンモデルのサポート
* exp2synシナプスモデルのサポート
* 形態ファイル（SWC）の自動前処理
  - Perisomatic modelへの変換
  - ゼロベースインデックスへの変換
  - 深さ優先探索（DFS）によるソート
* イオンチャネル設定の自動変換（JSON→CSV）
* MPI並列実行のサポート

**Neulite**

高性能シミュレータ：

* C言語で実装された軽量カーネル
* Allen Cell Types Databaseの Perisomatic modelに特化
* 定常電流入力のサポート
* ラズベリーパイから富岳まで幅広い環境で動作

**生成ファイル**

* `SONATA <https://github.com/AllenInstitute/sonata/tree/master>`_ 形式（BMTK互換）
* Neulite形式
  - ``<network_name>_population.csv`` - ニューロンポピュレーション定義
  - ``<src>_<trg>_connection.csv`` - 詳細なシナプス接続情報
  - 処理済みSWCファイル
  - イオンチャネル設定CSV
  - ``config.h`` - シミュレーション設定ヘッダー

**ドキュメント**

* プロジェクト概要
* bionet_liteとNeuliteの詳細説明
* セットアップガイド
* 4つのチュートリアル
  - Tutorial 01: 単一セルシミュレーション
  - Tutorial 02: スパイク入力の理解
  - Tutorial 03: 単一ポピュレーション
  - Tutorial 04: 複数ポピュレーション
* ユーザーガイド
  - 基本的な使い方
  - 設定ファイル
  - 並列実行
* アーキテクチャ
  - システム概要
  - 設計と実装（BMTKとBioNetの背景を含む）
* APIリファレンス
  - NeuliteBuilder API
  - ファイル形式仕様
* 高度なトピック
  - Allen V1モデルの実例
  - 仕様と制限事項
* FAQ

対応環境
--------

* Linux（推奨）
* macOS
* Python 3.7以上
* MPI環境（並列実行時）

既知の制限
----------

**対応モデル:**

* biophysicalモデルのみ対応（point neuron、virtual neuronは非対応）
* Perisomatic modelのみ対応

**対応シナプス:**

* exp2synのみ対応

**入力方式:**

* 定常電流入力のみ対応
* スパイク入力、Virtual cellによる入力は非対応

**その他:**

* delay値はint型のみ対応（小数点以下は丸められる）
* シナプス接続位置はネットワーク構築時にランダムに決定（bionetとは異なる）

今後の予定
==========

* スパイク入力機能の追加（Neulite側での実装）
* 新しいシナプスモデルのサポート
* パフォーマンスの継続的な改善
* ドキュメントの拡充
* チュートリアルの追加
