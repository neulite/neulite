#!/usr/bin/ruby

require 'fileutils'

files = Dir.glob("*.{c,h}") - ["config.h"]
p files

files.each{|file|
  FileUtils.cp(file, "#{file}.bak")
  open(file, "w"){|o|
    o.puts "// SPDX-License-Identifier: GPL-2.0-only"
    o.puts "// Copyright (C) 2024,2025 Neulite Core Team <neulite-core@numericalbrain.org>"
    o.puts ""
    IO.foreach("#{file}.bak"){|l| o.puts l }
  }
}

=begin
files = Dir.glob("*.bak")
p files
files.each{|file|
  orig_name = file.gsub(/.bak/,'')
  p orig_name
  FileUtils.mv(file, orig_name);
}
=end
