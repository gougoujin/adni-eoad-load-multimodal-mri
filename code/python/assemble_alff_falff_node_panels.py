"""
Repository usage summary
Description: Assemble ALFF and fALFF node-rendering images into a labeled two-panel composite figure.
Usage: Run from the command line with left image, right image, and output path arguments.
Outputs: Composite TIFF/PNG image files for manuscript or supplementary presentation.
Note: Update local input/output paths or pass explicit command-line arguments before running.

Public repository version: file paths are intentionally left blank or configurable.
This script was lightly sanitized for sharing and may require project-specific inputs.
"""

from PIL import Image, ImageDraw, ImageFont, ImageChops


def trim_white(img, threshold=248, pad=4):
    if img.mode != "RGB":
        img = img.convert("RGB")
    bg = Image.new("RGB", img.size, (255, 255, 255))
    diff = ImageChops.difference(img, bg).convert("L")
    mask = diff.point(lambda p: 255 if p > (255 - threshold) else 0)
    bbox = mask.getbbox()
    if bbox is None:
        return img
    l, t, r, b = bbox
    l = max(0, l - pad)
    t = max(0, t - pad)
    r = min(img.width, r + pad)
    b = min(img.height, b + pad)
    return img.crop((l, t, r, b))


def load_font(size=72):
    candidates = [
        "arialbd.ttf",
        "Arial Bold.ttf",
        "arial.ttf",
        "DejaVuSans-Bold.ttf",
    ]
    for fp in candidates:
        try:
            return ImageFont.truetype(fp, size=size)
        except:
            pass
    return ImageFont.load_default()


def assemble_alff_falff_node_panels(left_file=None, right_file=None, out_file=None):
    if not left_file or not right_file or not out_file:
        raise ValueError("left_file, right_file, and out_file must be provided.")

    A = Image.open(left_file).convert("RGB")
    B = Image.open(right_file).convert("RGB")

    A = trim_white(A, threshold=248, pad=4)
    B = trim_white(B, threshold=248, pad=4)

    target_h = max(A.height, B.height)
    if A.height != target_h:
        A = A.resize((round(A.width * target_h / A.height), target_h), Image.LANCZOS)
    if B.height != target_h:
        B = B.resize((round(B.width * target_h / B.height), target_h), Image.LANCZOS)

    gap = 36
    pad_left = 24
    pad_right = 24
    pad_bottom = 24
    title_band = 90

    canvas_w = pad_left + A.width + gap + B.width + pad_right
    canvas_h = title_band + target_h + pad_bottom
    canvas = Image.new("RGB", (canvas_w, canvas_h), (255, 255, 255))

    ax = pad_left
    ay = title_band
    bx = ax + A.width + gap
    by = title_band

    canvas.paste(A, (ax, ay))
    canvas.paste(B, (bx, by))

    draw = ImageDraw.Draw(canvas)
    font_panel = load_font(78)
    font_title = load_font(52)

    draw.text((ax + 0, 8), "A", font=font_panel, fill=(0, 0, 0))
    draw.text((bx + 0, 8), "B", font=font_panel, fill=(0, 0, 0))

    title_a = "ALFF significant nodes"
    title_b = "fALFF significant nodes"

    bbox_a = draw.textbbox((0, 0), title_a, font=font_title)
    bbox_b = draw.textbbox((0, 0), title_b, font=font_title)

    title_a_w = bbox_a[2] - bbox_a[0]
    title_b_w = bbox_b[2] - bbox_b[0]

    center_a_x = ax + A.width // 2
    center_b_x = bx + B.width // 2

    title_y = 28

    draw.text((center_a_x - title_a_w // 2, title_y), title_a, font=font_title, fill=(0, 0, 0))
    draw.text((center_b_x - title_b_w // 2, title_y), title_b, font=font_title, fill=(0, 0, 0))

    canvas.save(out_file, compression=None)
    out_png = out_file.rsplit(".", 1)[0] + ".png"
    canvas.save(out_png)

    print("Saved:")
    print(out_file)
    print(out_png)


def main():
    import argparse
    parser = argparse.ArgumentParser(description="Assemble ALFF and fALFF node-rendering images into a labeled composite panel.")
    parser.add_argument("--left-image", required=True, help="Path to the ALFF image.")
    parser.add_argument("--right-image", required=True, help="Path to the fALFF image.")
    parser.add_argument("--output", required=True, help="Output image path.")
    args = parser.parse_args()
    assemble_alff_falff_node_panels(args.left_image, args.right_image, args.output)

if __name__ == "__main__":
    main()
