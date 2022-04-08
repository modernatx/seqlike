from types import *
import pytest
import pandas as pd
from Bio.Align import MultipleSeqAlignment
from weblogo.color import Color
from seqlike.draw_utils import draw_alignment, view_alignment
from seqlike.draw_utils import convert_weblogo_color, apply_matching_colorscheme
from seqlike.draw_utils import aa_chemistry_simple, convert_colorscheme_to_color_map
from seqlike.SeqLike import SeqLike
from .test_SeqLike_alignment import get_aligned_aa_seqrecs_duplicate_ids


def test_convert_weblogo_color():

    color_format = "rgb"
    assert convert_weblogo_color(Color.from_string("red"), color_format) == (255, 0, 0)
    assert convert_weblogo_color(Color.from_string("green"), color_format) == (0, 128, 0)
    assert convert_weblogo_color(Color.from_string("blue"), color_format) == (0, 0, 255)
    assert convert_weblogo_color(Color.from_string("black"), color_format) == (0, 0, 0)
    assert convert_weblogo_color(Color.from_string("white"), color_format) == (255, 255, 255)

    color_format = "hex"
    assert convert_weblogo_color(Color.from_string("red"), color_format) == "#FF0000"
    assert convert_weblogo_color(Color.from_string("green"), color_format) == "#008000"
    assert convert_weblogo_color(Color.from_string("blue"), color_format) == "#0000FF"
    assert convert_weblogo_color(Color.from_string("black"), color_format) == "#000000"
    assert convert_weblogo_color(Color.from_string("white"), color_format) == "#FFFFFF"


def test_convert_colorscheme_to_color_map():
    # RGB
    color_dict = convert_colorscheme_to_color_map(aa_chemistry_simple, color_format="rgb")
    for letter in "GSTYCNQ":
        assert color_dict[letter] == (0, 128, 0)
    for letter in "KRH":
        assert color_dict[letter] == (0, 0, 255)
    for letter in "DE":
        assert color_dict[letter] == (255, 0, 0)
    for letter in "PAWFLIMV":
        assert color_dict[letter] == (0, 0, 0)
    for rule in aa_chemistry_simple().rules:
        for letter in rule.symbols:
            color_dict[letter] == convert_weblogo_color(rule.color, "rgb")
    # hex
    color_dict = convert_colorscheme_to_color_map(aa_chemistry_simple, color_format="hex")
    for rule in aa_chemistry_simple().rules:
        for letter in rule.symbols:
            color_dict[letter] == convert_weblogo_color(rule.color, color_format="hex")


def test_apply_matching_colorscheme():
    # gap_match
    assert apply_matching_colorscheme("-", "-", color_format="rgb") == (173, 216, 230)
    assert apply_matching_colorscheme("-", "-", color_format="hex") == "#ADD8E6"
    # gap
    assert apply_matching_colorscheme("-", "S", color_format="rgb") == (255, 255, 255)
    assert apply_matching_colorscheme("-", "S", color_format="hex") == "#FFFFFF"
    # match
    assert apply_matching_colorscheme("S", "S", color_format="rgb") == (50, 205, 50)
    assert apply_matching_colorscheme("S", "S", color_format="hex") == "#32CD32"
    # mismatch
    assert apply_matching_colorscheme("S", "A", color_format="rgb") == (139, 0, 0)
    assert apply_matching_colorscheme("S", "A", color_format="hex") == "#8B0000"


def test_draw_alignment():
    aligned = get_aligned_aa_seqrecs_duplicate_ids()
    fig = draw_alignment(MultipleSeqAlignment(aligned))


def test_view_alignment():
    try:
        import bokeh as bk
    except ImportError:
        pytest.skip("Bokeh not installed, skipping view_alignment() test")
    try:
        import IPython.display
    except ImportError:
        pytest.skip("IPython not installed, skipping view_alignment() test")
    aligned = get_aligned_aa_seqrecs_duplicate_ids()
    fig = view_alignment(MultipleSeqAlignment(aligned))


@pytest.mark.xfail(reason="Will fail is Ghostscript is not on path.")
def test_draw_weblogo():
    aligned = get_aligned_aa_seqrecs_duplicate_ids()
    df = pd.DataFrame({"seqs": [SeqLike(seq, seq_type="aa") for seq in aligned]})
    fig = df.seqs.seq.weblogo()
