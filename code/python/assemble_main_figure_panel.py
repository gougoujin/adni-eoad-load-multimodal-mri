"""
Repository usage summary
Description: Assemble the main manuscript panel from a left brain image and a right summary panel image.
Usage: Run from the command line with explicit input image paths and an output path.
Outputs: Composite PNG/TIFF image files for the final assembled panel.
Note: Update local input/output paths or pass explicit command-line arguments before running.

Public repository version: file paths are intentionally left blank or configurable.
This script was lightly sanitized for sharing and may require project-specific inputs.
"""

from PIL import Image, ImageOps, ImageChops, ImageDraw, ImageFont


def trim_white(img, bg_threshold=248, pad=4):
    if img.mode != "RGB":
        img = img.convert("RGB")
    bg = Image.new("RGB", img.size, (255, 255, 255))
    diff = ImageChops.difference(img, bg).convert("L")
    bbox = diff.point(lambda p: 255 if p > (255 - bg_threshold) else 0).getbbox()
    if bbox is None:
        return img
    left, top, right, bottom = bbox
    left = max(0, left - pad)
    top = max(0, top - pad)
    right = min(img.width, right + pad)
    bottom = min(img.height, bottom + pad)
    return img.crop((left, top, right, bottom))


def load_font(size=110):
    font_candidates = [
        "arialbd.ttf",
        "Arial Bold.ttf",
        "arial.ttf",
        "DejaVuSans-Bold.ttf",
    ]
    for fp in font_candidates:
        try:
            return ImageFont.truetype(fp, size=size)
        except:
            pass
    return ImageFont.load_default()


def assemble_fig1_main_v3(brain_img_file=None, panel_b_img_file=None, out_file=None):
    if not brain_img_file or not panel_b_img_file or not out_file:
        raise ValueError("brain_img_file, panel_b_img_file, and out_file must be provided.")

    A = Image.open(brain_img_file).convert("RGB")
    B = Image.open(panel_b_img_file).convert("RGB")

    # 裁掉外围空白，避免 A 因为留白太多显得变小
    A = trim_white(A, bg_threshold=248, pad=4)
    B = trim_white(B, bg_threshold=248, pad=4)

    # 统一高度，但保持各自比例
    target_h = max(A.height, B.height)
    if A.height != target_h:
        A = A.resize((round(A.width * target_h / A.height), target_h), Image.LANCZOS)
    if B.height != target_h:
        B = B.resize((round(B.width * target_h / B.height), target_h), Image.LANCZOS)

    # 左边额外扩张，专门给 A 放位置
    left_expand = 130
    gap = 28
    pad_top = 18
    pad_bottom = 18
    pad_right = 24

    canvas_w = left_expand + A.width + gap + B.width + pad_right
    canvas_h = target_h + pad_top + pad_bottom
    canvas = Image.new("RGB", (canvas_w, canvas_h), (255, 255, 255))

    # 垂直居中
    a_x = left_expand
    a_y = (canvas_h - A.height) // 2
    b_x = a_x + A.width + gap
    b_y = (canvas_h - B.height) // 2

    canvas.paste(A, (a_x, a_y))
    canvas.paste(B, (b_x, b_y))

    # 画 A，字号明显调大，尽量和 B/C/D/E 接近
    draw = ImageDraw.Draw(canvas)
    font = load_font(size=110)
    text_x = 18
    text_y = max(6, a_y + 2)
    draw.text((text_x, text_y), "A", font=font, fill=(0, 0, 0))

    # 保存 png 和 tif，同名覆盖
    canvas.save(out_file)
    tif_file = out_file.rsplit(".", 1)[0] + ".tif"
    canvas.save(tif_file, compression=None)

    print("Saved:")
    print(out_file)
    print(tif_file)


def main():
    import argparse
    parser = argparse.ArgumentParser(description="Assemble a two-panel main figure from existing image files.")
    parser.add_argument("--brain-image", required=True, help="Path to the left image file.")
    parser.add_argument("--panel-b-image", required=True, help="Path to the right image file.")
    parser.add_argument("--output", required=True, help="Output image path.")
    args = parser.parse_args()
    assemble_fig1_main_v3(args.brain_image, args.panel_b_image, args.output)

if __name__ == "__main__":
    main()
