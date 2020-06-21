Sphinx 样例
##############################


简介
----

学习一下Sphinx，制作一个PDF版本。

相关
----

- 书的文字部分来自于 `Andriki`_ 提供的Mediawiki源码；
- 使用 `Sphinx`_ 制作文档
- 配置`【杂谈】windows10配置make命令`_
- 了解到Sphinx这一文档写作工具，一个案例`跟我一起写Makefile重制版`_

本地编译
--------

#. Clone项目到本地::

   $ git clone https://git.dev.tencent.com/lingr7/Sphinx_quickstart_test.git

#. 安装依赖::

   $ pip install -r requirements.txt

#. 编译生成HTML::

   $ make html
   $ firefox build/html/index.html&

#. 编译生成PDF（要求安装TeXLive 2016）::

   $ make latexpdf
   $ evince build/latex/quasi-newtonmethod.pdf&

.. _`陈皓`: http://coolshell.cn/haoel
.. _`Andriki`: http://andriki.com/mediawiki/index.php?title=Linux:%E8%B7%9F%E6%88%91%E4%B8%80%E8%B5%B7%E5%86%99Makefile
.. _`Sphinx`: http://sphinx-doc.org/
.. _`GNU Make Manual`: https://www.gnu.org/software/make/manual/
.. _`跟我一起写Makefile重制版`:https://github.com/seisman/how-to-write-makefile
.. _`【杂谈】windows10配置make命令`:https://blog.csdn.net/C2681595858/article/details/85554359