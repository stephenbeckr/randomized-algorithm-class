# Homeworks

For APPM 5650 Randomized Algorithms, Fall 2021.

HW solutions are on Canvas

- [HW 1](APPM5650Fall21_RandomizedAlgos_HW01.pdf), due Friday Aug 27, 2021 (topics: read some of Mahoney, Martinsson and Tropp, and Vershynin)
- [HW 2](APPM5650Fall21_RandomizedAlgos_HW02.pdf), due Monday Sep 6, 2021 (topics: linear algebra, sparse matrices, Freivald's algorithm, random orthogonal matrices). Turn this in via Gradescope since there's no class Monday (Labor Day)
- [HW 3](APPM5650Fall21_RandomizedAlgos_HW03.pdf), due Monday Sep 13, 2021 (topics: linear algebra, basis probability, account on [CU's research computing](https://www.colorado.edu/rc).
- [HW 4](APPM5650Fall21_RandomizedAlgos_HW04.pdf), due Monday Sep 20, 2021 (topics: large data files, random projections vs tSNE)
- [HW 5](APPM5650Fall21_RandomizedAlgos_HW05.pdf), due Monday Sep 27, 2021 (topics: random projections)
- [HW 6](APPM5650Fall21_RandomizedAlgos_HW06.pdf), due Monday Oct 4, 2021 (topics: different types of sampling without replacement)

# Turning in HW
While the class meets in person, please turn in **hard-copies** of your homework in class.

If we end up meeting remotely due to COVID, then we will switch to gradescope.  As of now, we are *not using Gradescope*.


## FAQ for electronic submissions

### General

**Gradescope** has a [submission guide](https://gradescope-static-assets.s3.amazonaws.com/help/submitting_hw_guide.pdf) that recommends software for your phone to take pictures of written homework and convert it to a PDF (your final submission to Gradescope must be a PDF).

Note: the links in the PDFs will not work if you view the PDF on github, but if you open the PDF in its own tab, or download it, all the links should work.

**Collaboration**: Collaboration with your fellow students is OK and in fact recommended, although direct copying is not allowed.  Please write down the names of the students that you worked with.

**Internet**: The internet is allowed for basic tasks (e.g., looking up definitions on wikipedia) but it is
not permissible to search for proofs or to post requests for help on forums such as [math.stackexchange.com](http://math.stackexchange.com/)
or to look at solution manuals

#### Merging multiple PDF files

**Mac** You can use the Preview software that comes with Mac, and drag-and-drop in the Thumbnail view, or follow these [instructions](https://support.apple.com/en-us/HT202945).

**Linux** install `pdftk` (e.g., `apt-get install pdftk`), and the on the command line, it's just `pdftk inputFile1.pdf inputFile2.pdf cat output outputFileName.pdf`.  This works on Mac and Windows too (on Mac, the exact command line works; on Windows, I'm not sure).

**Windows** there are [lists of free web- and desktop-based software](https://superuser.com/a/34294), but [PDFtk](https://www.pdflabs.com/tools/pdftk-the-pdf-toolkit/) is one of the most classic and respected (no viruses). I haven't used PDFtk on Windows, but the website claims they have a GUI; or if you don't like their GUI, try a [3rd party GUI that uses PDFtk](https://www.pdflabs.com/tools/pdftk-the-pdf-toolkit/).

### Python
For overall Python, and numpy in particular, Matlab users might like [NumPy for Matlab users](https://numpy.org/doc/stable/user/numpy-for-matlab-users.html).

For **plotting** in Python using Matplotlib, try these [plotting cheatsheets](https://github.com/matplotlib/cheatsheets) and [controlling figure aesthetics](https://seaborn.pydata.org/tutorial/aesthetics.html) with seaborn.


#### Jupyter

Tips for exporting jupyter notebook code to a PDF:

- You can try this [Notebook to PDF conversion website](https://htmtopdf.herokuapp.com/ipynbviewer/) that some of our students have had good luck with

- Or try `nbconvert` which requires [`pandoc`](https://pandoc.org/installing.html). You can do this on [Colab](https://colab.research.google.com/), following the [instructions here](https://stackoverflow.com/a/54191922) (but note that you may need to add a backslash before any white space when you run commands, e.g., change a command like

`!cp drive/My Drive/Colab Notebooks/Untitled.ipynb ./`
to
``!cp drive/My\ Drive/Colab\ Notebooks/Untitled.ipynb ./``
)

Note that if you include latex in the jupyter notebook, when you run `nbconvert`, you cannot have any whitespace near the `\$` symbols for math due to a requirement of `pandoc` (see [here](https://pandoc.org/MANUAL.html#extension-tex_math_dollars)).  So, ``$ f(x) = 3x^2 $`` will not work, but `$f(x) = 3x^2$` will be OK.

The downside of `nbconvert` is that images are saved as png, not pdf, so fonts don't come through, but that's not a big deal for homework.

If you run jupyter locally, you might be able to run `nbconvert` without using the command line; go to "Download" the "PDF via LaTeX".

#### Python source code (not Jupyter)
The *non-preferred* ways are (1) screenshot of your editor (not so nice since it's an image not text, but at least you get syntax color highlighting), and (2) export from a text editor to PDF (not so nice if you don't get syntax color highlighting).

It's not nice to the graders to submit code without syntax color highlighting!

Better ways: it depends on your system and editor, but there are many ways. For example, this [stackoverflow 'printing python code to PDF'](https://stackoverflow.com/q/20412038) offers several suggestions. For example, since I already use `vim` and its setup with syntax highlighting, I can do [this answer](https://stackoverflow.com/a/20412421) and do `vim abc.py -c ":hardcopy > abc.ps" -c ":q"` followed by `ps2pdf abc.ps abc.pdf` -- no extra software needed!


### Matlab

You can use the export notebook features in Matlab (it can handle latex) if you want; see the [Live-Editor](https://www.mathworks.com/products/matlab/live-editor.html); there are also claims on the internet that it's easy to get Jupyter to run with a Matlab kernel, so you could use Jupyter.

To just export a figure, there are builtin methods, but one of the nicer ways is to use [export_fig](https://www.mathworks.com/matlabcentral/fileexchange/23629-export_fig), which works like `export_fig MyFileName -pdf -transparent` and makes a file `MyFileName.pdf` (note that PDF files for figures are preferred, since then the text is saved as a font and not bitmapped)
