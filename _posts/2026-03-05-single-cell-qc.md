---
title: "单细胞测序数据质控"
date: 2026-03-05
permalink: /posts/2026/03/单细胞测序数据质控/
tags:
  - Bioinformatics
  - Single Cell
  - QC
---

# https://mp.weixin.qq.com/s/paEADUSr5dzHXNYQsSjz6g

![cover_image](https://mmbiz.qpic.cn/mmbiz_jpg/zd9Yq54C2EESKvMTWbicEXlI3SVPLCDV7UZ6KQlvg2dr904DcmBd9vibsl3fqI8FNgAOMic8aOiavkSDAhOlTSAa8g/0?wx_fmt=jpeg)

# Rust学习笔记-01

![image]()

## 前言-为什么想学Rust？

一直想多学一门编译语言，主要有几个原因吧（1. 看到一位老师实验室要求需要掌握一门编译语言；2. 自己想试着开发一些实用的生信工具，现阶段自己现在相比于数据分析，探索生物学层面的意义而言更想做一些实用性的东西出来），之前在C、C++、Rust之间摇摆，一直没有行动起来学，确实也怪自己拖延。

直到前段时间，想去的那个实验室发了一篇很有创新性的文章，又勾起了我想前往深造的想法。（叠个甲，不是只想着为了发文章才去，而是觉得在那里会有更多的可能性。）在浏览实验室信息的时候，偶然发现一位博士师兄竟然是我之前在B-G-I实习时候的室友，尽管只当了三四天的室友，但对他的一些基本情况还是比较了解的。后来知道他去了那个学校，但没想到竟然去了那位老师门下。

又在网上搜索这位师兄的相关信息，确实很优秀，掌握很多技能，渐渐地心里将他变成了我研究生阶段的榜样。看到他的github或者是推文中都在推荐Rust，我那摇摆的天平最终偏向Rust这边，下决心要开始行动了。

再不做真的就wan（晚/完）了。

说干就干，赶紧学！！！

PS：现在脑中有个小工具的雏形，争取能在今年实现！！

## Rust资源

https://course.rs/about-book.html

https://www.runoob.com/rust/rust-tutorial.html

https://rustwiki.org/

https://github.com/rust-lang/rustlings/

https://play.rust-lang.org/?version=stable&mode=debug&edition=2024

https://doc.rust-lang.org/stable/rust-by-example/index.html

## Rust特点

简单了解一下Rust特点，之前一直听说一大特点是内存安全，之前在大一时候学过C++，但是到后面学指针的时候就开小差了，其实也不能很好理解这部分内容。这里引用一下菜鸟教程里的内容吧：

对于生信可以借助Rust进行工具开发，这也是我一直想学着做的一件事

## Cargo

我们在平常管理生信环境是经常会用到conda/pip，那么在Rust里面可以使用Cargo来对包进行管理。

### 创建项目

在这里使用cargo创建一个world\_hello项目

`(base) PS D:\000zyf\Learning\rust_learn> cargo new world_hello  
    Creating binary (application) `world_hello` package  
note: see more `Cargo.toml` keys and their definitions at https://doc.rust-lang.org/cargo/reference/manifest.html  
(base) PS D:\000zyf\Learning\rust_learn> cd .\world_hello\  
(base) PS D:\000zyf\Learning\rust_learn\world_hello> ls  
    Directory: D:\000zyf\Learning\rust_learn\world_hello  
Mode                 LastWriteTime         Length Name  
----                 -------------         ------ ----  
d----           2025/4/12    21:33                src  
-a---           2025/4/12    21:33             82 Cargo.toml  
(base) PS D:\000zyf\Learning\rust_learn\world_hello>`

这里如果使用命令cargo new --vcs git xxx来创建项目的话还会给出相应的git文件。

吐槽：windows终端上操作真难受！！

### 运行项目

有两种方式运行：

#### cargo run

在刚刚创建的项目路径下运行cargo run便可看到：

`(base) PS D:\000zyf\Learning\rust_learn\world_hello> cargo run  
  
   Compiling world_hello v0.1.0 (D:\000zyf\Learning\rust_learn\world_hello)  
    Finished `dev` profile [unoptimized + debuginfo] target(s) in 5.83s  
     Running `target\debug\world_hello.exe`  
Hello, world!`

上述代码，cargo run 首先对项目进行编译，然后再运行，因此它实际上等同于运行了两个指令

#### 手动编译运行

`(base) PS D:\000zyf\Learning\rust_learn\world_hello> cargo build  
    Finished `dev` profile [unoptimized + debuginfo] target(s) in 0.03s  
(base) PS D:\000zyf\Learning\rust_learn\world_hello> ls  
    Directory: D:\000zyf\Learning\rust_learn\world_hello  
Mode                 LastWriteTime         Length Name    
----                 -------------         ------ ----    
d----           2025/4/12    21:39                src         
d----           2025/4/12    21:50                target      
-a---           2025/4/12    21:39              8 .gitignore  
-a---           2025/4/12    21:50            155 Cargo.lock  
-a---           2025/4/12    21:39             82 Cargo.toml  
(base) PS D:\000zyf\Learning\rust_learn\world_hello> .\target\debug\world_hello.exe  
Hello, world!  
(base) PS D:\000zyf\Learning\rust_learn\world_hello>`

行云流水，但谈不上一气呵成。 细心的读者可能已经发现，在调用的时候，路径 ./target/debug/world\_hello 中有一个明晃晃的 debug 字段，没错我们运行的是 debug 模式，在这种模式下，代码的编译速度会非常快，可是福兮祸所伏，运行速度就慢了. 原因是，在 debug 模式下，Rust 编译器不会做任何的优化，只为了尽快的编译完成，让你的开发流程更加顺畅。

比如你想要高性能的代码怎么办？ 简单，添加 --release 来编译：

cargo run --release

cargo build --release

试着运行一下我们高性能的 release 程序：

`(base) PS D:\000zyf\Learning\rust_learn\world_hello> cd .\target\           
(base) PS D:\000zyf\Learning\rust_learn\world_hello\target> ls  
    Directory: D:\000zyf\Learning\rust_learn\world_hello\target  
Mode                 LastWriteTime         Length Name  
----                 -------------         ------ ----  
d----           2025/4/12    21:51                debug  
d----           2025/4/12    22:00                release  
-a---           2025/4/12    21:51           1120 .rustc_info.json  
-a---           2025/4/12    21:50            177 CACHEDIR.TAG  
(base) PS D:\000zyf\Learning\rust_learn\world_hello\target> .\release\world_hello.exe  
Hello, world!`

### cargo check

当项目大了后，cargo run 和 cargo build 不可避免的会变慢，那么有没有更快的方式来验证代码的正确性呢？大杀器来了，接着！

cargo check 是我们在代码开发过程中最常用的命令，它的作用很简单：快速的检查一下代码能否编译通过。因此该命令速度会非常快，能节省大量的编译时间。

### Cargo.toml和Cargo.lock

Cargo.toml 和 Cargo.lock 是 cargo 的核心文件，它的所有活动均基于此二者。

Cargo.toml 是 cargo 特有的项目数据描述文件。它存储了项目的所有元配置信息，如果 Rust 开发者希望 Rust 项目能够按照期望的方式进行构建、测试和运行，那么，必须按照合理的方式构建 Cargo.toml。

Cargo.lock 文件是 cargo 工具根据同一项目的 toml 文件生成的项目依赖详细清单，因此我们一般不用修改它。

什么情况下该把 Cargo.lock 上传到 git 仓库里？很简单，当你的项目是一个可运行的程序时，就上传 Cargo.lock，如果是一个依赖库项目，那么请把它添加到 .gitignore 中。

现在用 VSCode 打开上面创建的"世界，你好"项目，然后进入根目录的 Cargo.toml 文件，可以看到该文件包含不少信息：

`[package]  
name = "world_hello"  
version = "0.1.0"  
edition = "2024"  
[dependencies]`

name 字段定义了项目名称，version 字段定义当前版本，新项目默认是 0.1.0，edition 字段定义了我们使用的 Rust 大版本。

### 定义项目依赖

使用 cargo 工具的最大优势就在于，能够对该项目的各种依赖项进行方便、统一和灵活的管理。

在 Cargo.toml 中，主要通过各种依赖段落来描述该项目的各种依赖项：

# Rust语言圣经(Rust Course)：https://course.rs/first-try/cargo.html

![image]()

微信扫一扫  
关注该公众号

![](http://mmbiz.qpic.cn/mmbiz_png/zd9Yq54C2EHDEW8BRcCexuLYsvra6xFaIBZIyS46JXuDibpiaBDFxrqYXp4u09HPmFhC5jqcS9J8M3CPiaykRNgkA/0?wx_fmt=png)
![image]()
![跳转二维码]()
![作者头像](http://mmbiz.qpic.cn/mmbiz_png/zd9Yq54C2EHDEW8BRcCexuLYsvra6xFaIBZIyS46JXuDibpiaBDFxrqYXp4u09HPmFhC5jqcS9J8M3CPiaykRNgkA/0?wx_fmt=png)

微信扫一扫可打开此内容，  
使用完整服务

---

