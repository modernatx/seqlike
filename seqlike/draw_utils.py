import sys
import os
import numpy as np
from typing import Callable, Union

from PIL import Image, ImageDraw, ImageFont

import lazy_loader as lazy

wl = lazy.load("weblogo")

# try:
# import bokeh as bk
bk = lazy.load("bokeh")
# from bokeh.plotting import figure, show
# from bokeh.core.properties import value

# for visualization in jupyter notebook
try:
    get_ipython
    from bokeh.io import output_notebook

    output_notebook()
except NameError:
    pass
except ImportError:
    raise ValueError(
        "Bokeh is required for the interactive alignment viewer in notebooks. Install `seqlike[notebook]`",
    )

from .alphabets import gap_letter


# revised from chemistry_extended: removed neutral and moved N and Q to polar
def aa_chemistry_simple():
    return wl.colorscheme.ColorScheme(
        [
            wl.colorscheme.SymbolColor("GSTYCNQ", "green", "polar"),
            wl.colorscheme.SymbolColor("KRH", "blue", "basic"),
            wl.colorscheme.SymbolColor("DE", "red", "acidic"),
            wl.colorscheme.SymbolColor("PAWFLIMV", "black", "hydrophobic"),
            wl.colorscheme.SymbolColor("X", "gray", "unknown"),
        ],
        alphabet=wl.seq.protein_alphabet,
    )


# from makelogo
def aa_chemistry_extended():
    return wl.colorscheme.ColorScheme(
        [
            wl.colorscheme.SymbolColor("GSTYC", "green", "polar"),
            wl.colorscheme.SymbolColor("NQ", "purple", "neutral"),
            wl.colorscheme.SymbolColor("KRH", "blue", "basic"),
            wl.colorscheme.SymbolColor("DE", "red", "acidic"),
            wl.colorscheme.SymbolColor("PAWFLIMV", "black", "hydrophobic"),
            wl.colorscheme.SymbolColor("X", "gray", "unknown"),
        ],
        alphabet=wl.seq.protein_alphabet,
    )


def aa_dssp_color():
    return wl.colorscheme.ColorScheme(
        [
            wl.colorscheme.SymbolColor("EB", "red", "strand"),
            wl.colorscheme.SymbolColor("HGI", "blue", "helix"),
            wl.colorscheme.SymbolColor("TSC-", "light gray", "coil"),
        ],
        alphabet=wl.seq.protein_alphabet,
    )


def nt_simple():
    return wl.colorscheme.ColorScheme(
        [
            wl.colorscheme.SymbolColor("G", "orange"),
            wl.colorscheme.SymbolColor("TU", "red"),
            wl.colorscheme.SymbolColor("C", "blue"),
            wl.colorscheme.SymbolColor("A", "green"),
        ],
    )


def convert_weblogo_color(color: "wl.color.Color", color_format: str) -> Union[tuple, str]:
    """Convert weblogo Color to Bokeh color object

    Note: Weblogo colors are RGB but fractional [0, 1],
    whereas Bokeh and draw_alignment are [0, 255]

    :sa: https://github.com/WebLogo/weblogo/blob/master/weblogo/color.py

    :param color: A weblogo Color object.
    :param color_format: Either "rgb" or "hex".
    :returns: Either an RGB tuple (for "rgb") or hexadecimal string (for "hex").
    """
    assert color_format in ["rgb", "hex"]

    rgb_tuple = int(255 * color.red), int(255 * color.green), int(255 * color.blue)
    hex_str = f"#{rgb_tuple[0]:02X}{rgb_tuple[1]:02X}{rgb_tuple[2]:02X}"

    if color_format == "rgb":
        return rgb_tuple
    else:
        return hex_str


def convert_colorscheme_to_color_map(color_scheme: Callable, color_format: str) -> dict:
    """Convert weblogo ColorScheme into bokeh color map
    :param color_scheme: a Callable that returns a weblogo ColorScheme object
    :param color_format: 'hex' or 'rgb' for hex string or RGB tuple, respectively
    :returns: a dict of bokeh colors indexed by letter
    """
    # return a Bokeh color object or simple RGB tuple (for draw_alignment)

    # convert SymbolColor to bokeh color object
    color_dict = dict()
    for rule in color_scheme().rules:
        color_dict[rule.symbols] = convert_weblogo_color(rule.color, color_format)
    # default for spaces (white)
    color_dict["-*"] = convert_weblogo_color(wl.color.Color.from_string("white"), color_format)
    # expand letter strings so that dict maps to single letters
    expanded_color_dict = dict()
    for letters, color in color_dict.items():
        expanded_color_dict.update(dict((l, color) for l in letters))
    return expanded_color_dict


def apply_matching_colorscheme(letter, ref_letter, color_format: str):
    """Apply a match/mismatch color scheme to sequence letters
    :param letter: letter from target sequence
    :param ref_letter: letter from reference sequence (for match/mismatch)
    :param color_format: 'hex' or 'rgb' for hex string or RGB tuple, respectively
    :returns: an RGB hex string (for Bokeh) or simple RGB tuple (for vizqes)
    """
    # gap match
    if letter == gap_letter and letter == ref_letter:
        return convert_weblogo_color(wl.color.Color.from_string("lightblue"), color_format)
    # gap
    elif letter == gap_letter:
        return convert_weblogo_color(wl.color.Color.from_string("white"), color_format)
    # match
    elif letter == ref_letter:
        return convert_weblogo_color(wl.color.Color.from_string("limegreen"), color_format)
    # mismatch
    else:
        return convert_weblogo_color(wl.color.Color.from_string("darkred"), color_format)


def find_font(size, fontpath=None):
    """Find and scale font based on fontpath.

    Helper function for draw_alignment.

    :param size: desired font size
    :param fontpath: optional search path for font file (.ttf)
    :returns: PIL.ImageFont object

    :sa: vizqes (https://pypi.python.org/pypi/vizqes)
    """
    if fontpath:
        font_searchpath = fontpath
    else:
        font_searchpath = os.path.join(os.path.dirname(__file__), "FreeMono.ttf")
    try:
        font = ImageFont.truetype(font_searchpath, size=size)
        sys.stdout.write("Found font in {}\n".format(str(font_searchpath)))
    except IOError as e:
        sys.stderr.write(str(e))
        sys.stderr.write("could not find font in {}\nUsing default\n".format(str(font_searchpath)))
        font = ImageFont.load_default()
    return font


def draw_alignment(
    aligned,
    colorscheme: Callable = aa_chemistry_simple,
    boxwidth=2,
    boxheight=12,
    label_width=100,
    show_ids=False,
    show_names=False,
    show_descriptions=False,
    show_grouping=False,
):
    """Generate a colored figure from an alignment
    :param aligned: MultipleSeqAlignment object
    :param colorscheme: a Callable that returns a weblogo ColorScheme object
    :param boxwidth: column width of alignment
    :param boxheight: row height of alignment
    :param label_width: maximum length of row label; if None, extend to maximum label length
    :param show_ids: if True, show SeqRecord ID for each row
    :param show_names: if True, show SeqRecord name for each row
    :param show_descriptions: if True, show SeqRecord description for each row
    :param show_grouping: if True, highlight changes from reference in red against green background,
        instead of using the residue colorscheme
    :returns: PIL Image object

    :note: based on vizqespkg.vizqes_main.draw
    :sa: vizqespkg.vizqes_main.draw
    :sa: http://www.bioinformatics.nl/~berndb/aacolour.html
    """

    if show_names or show_ids or show_descriptions:
        font = find_font(boxheight)
        offset = -1
        if show_names:
            offset += font.getsize(max([m.name[None:label_width] for m in aligned], key=len))[0] + 1
        if show_ids:
            offset += font.getsize(max([m.id[None:label_width] for m in aligned], key=len))[0] + 1
        if show_descriptions:
            offset += font.getsize(max([m.description[None:label_width] for m in aligned], key=len))[0] + 1
    else:
        font, offset = None, 0

    height = len(aligned) * boxheight
    width = aligned.get_alignment_length() * boxwidth + offset

    img = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(img)
    yd = None

    color_dict = convert_colorscheme_to_color_map(colorscheme, color_format="rgb")
    refseq = aligned[0].seq
    for y, member in enumerate(aligned):
        y *= boxheight
        for x, xs in enumerate(member.seq):
            if show_grouping:
                color = apply_matching_colorscheme(xs, refseq[x], color_format="rgb")
            else:
                color = color_dict[xs]
            x *= boxwidth
            for i in range(0, boxwidth):
                xd = x + i + offset
                for j in range(0, boxheight):
                    yd = y + j
                    draw.point((xd, yd), fill=color)
        if show_names or show_ids or show_descriptions:
            text = ""
            if show_names:
                text += member.name[None:label_width] + " "
            if show_ids:
                text += member.id[None:label_width] + " "
            if show_descriptions:
                text += member.description[None:label_width] + " "
            # clip last ' ' from text
            draw.text((0, yd - boxheight), text[:-1], font=font, fill=(0, 0, 0))
    return img


def view_alignment(
    aligned,
    fontsize="9pt",
    show_N=100,
    colorscheme: Callable = aa_chemistry_simple,
    boxwidth=9,
    boxheight=15,
    label_width=None,
    show_descriptions=False,
    show_grouping=False,
):
    """Bokeh sequence alignment view for protein and nucleic acid sequences


    :sa: https://dmnfarrell.github.io/bioinformatics/bokeh-sequence-aligner

    :param aligned: MultipleSeqAlignment object
    :param fontsize: font size for text labels
    :param show_N: size of sequence window (in number of sequence letters)
    :param colorscheme: a Callable that returns a weblogo ColorScheme object
    :param boxwidth: column width of alignment
    :param boxheight: row height of alignment
    :param label_width: maximum length of row label; if None, extend to maximum label length
    :param show_descriptions: if True, show SeqRecord description for each row
    :param show_grouping: if True, highlight changes from reference in red against green background,
        instead of using the residue colorscheme
    :returns: A Bokeh plot of the Multiple Sequence Alignment.
    """
    from bokeh.models import ColumnDataSource, Range1d
    from bokeh.plotting import figure
    from bokeh.models.glyphs import Rect

    def get_colors(seqs, color_scheme):
        """make colors for letters in sequence

        :param seqs: A string sequence.
        :param color_scheme: A string.
        :returns: a sequence of colors for each letter in seqs.
        """
        # get colors
        color_dict = convert_colorscheme_to_color_map(color_scheme, color_format="hex")
        # assign colors to sequences
        text = [i for s in list(seqs) for i in s]
        return [color_dict[a] for a in text]

    def get_colors_for_matching(seqs):
        """match/mismatch color scheme for show_grouping

        :param seqs: Sequences for which colors need to be matched.
        :returns: a list of colors (strings)
        """
        refseq = seqs[0]
        colors = list()
        for seq in list(seqs):
            for xs, ref_s in zip(seq, refseq):
                colors.append(apply_matching_colorscheme(xs, ref_s, color_format="hex"))
        return colors

    # make sequence and id lists from the aligned object
    seqs = [rec.seq for rec in (aligned)]
    if show_descriptions:
        labels = [f"{row} - {rec.description} ({rec.id})" for (row, rec) in enumerate(aligned)]
    else:
        labels = [f"{row} - {rec.id}" for (row, rec) in enumerate(aligned)]

    if label_width:
        labels = [label[:label_width] for label in labels]
    else:
        label_width = max(len(label) for label in labels)

    text = [i for s in list(seqs) for i in s]
    if show_grouping:
        colors = get_colors_for_matching(seqs)
    else:
        colors = get_colors(seqs, colorscheme)
    N = len(seqs[0])
    S = len(seqs)

    x = np.arange(1, N + 1)
    # need to reverse y so that sequences are plotted top-to-bottom
    y = np.arange(S - 1, -1, -1)
    # creates a 2D grid of coords from the 1D arrays
    xx, yy = np.meshgrid(x, y)
    # flattens the arrays
    gx = xx.ravel()
    gy = yy.flatten()
    # use recty for rect coords with an offset
    recty = gy + 0.5
    # now we can create the ColumnDataSource with all the arrays
    source = ColumnDataSource(dict(x=gx, y=gy, recty=recty, text=text, colors=colors))
    plot_height = len(seqs) * boxheight + 50
    x_range = Range1d(0, N + 1, bounds="auto")
    viewlen = min(show_N, N)
    # view_range is for the close up view
    view_range = (0, viewlen)
    tools = "xpan,xwheel_zoom,reset,save"

    # plot_width combines length of text labels and number of letters in sequence view window
    # note: this part requires additional tuning; 5 pixel average width of y-axis labels is a guess
    plot_width = int(5 * label_width) + boxwidth * viewlen + 40

    # entire sequence view (no text, with zoom)

    p = figure(
        title=None,
        width=plot_width,
        height=50,
        x_range=x_range,
        y_range=(0, S),
        tools=tools,
        min_border=0,
        toolbar_location="below",
    )
    rects = Rect(
        x="x",
        y="recty",
        width=1,
        height=1,
        fill_color="colors",
        line_color=None,
        fill_alpha=0.6,
    )
    p.add_glyph(source, rects)
    p.yaxis.visible = False
    p.grid.visible = False

    # sequence text view with ability to scroll along x axis
    p1 = figure(
        title=None,
        width=plot_width,
        height=plot_height,
        x_range=view_range,
        y_range=labels[::-1],
        tools="xpan,reset,save",
        min_border=0,
        toolbar_location="below",
    )  # , lod_factor=1)
    glyph = bk.models.glyphs.Text(
        x="x",
        y="y",
        text="text",
        text_align="center",
        text_color="black",
        text_font=bk.core.properties.value("monospace"),
        text_font_size=fontsize,
    )
    rects = Rect(
        x="x",
        y="recty",
        width=1,
        height=1,
        fill_color="colors",
        line_color=None,
        fill_alpha=0.4,
    )
    p1.add_glyph(source, glyph)
    p1.add_glyph(source, rects)

    p1.grid.visible = False
    p1.xaxis.major_label_text_font_style = "bold"
    p1.yaxis.minor_tick_line_width = 0
    p1.yaxis.major_tick_line_width = 0

    p = bk.layouts.gridplot([[p], [p1]], toolbar_location="below")
    bk.plotting.show(p)
    return p
