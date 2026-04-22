"""
Repository usage summary
Description: Stack two pre-rendered graph-summary images vertically into a single composite panel.
Usage: Run from the command line with top image, bottom image, and output path arguments.
Outputs: Composite PNG/TIFF image files containing the vertically stacked panels.
Note: Update local input/output paths or pass explicit command-line arguments before running.

Public repository version: file paths are intentionally left blank or configurable.
This script was lightly sanitized for sharing and may require project-specific inputs.
"""

from PIL import Image, ImageChops

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

def stack_graph_summary_panels(top_image=None, bottom_image=None, out_file=None):
    if not top_image or not bottom_image or not out_file:
        raise ValueError("top_image, bottom_image, and out_file must be provided.")

    fig2 = Image.open(top_image).convert("RGB")
    fig3 = Image.open(bottom_image).convert("RGB")

    fig2 = trim_white(fig2, threshold=248, pad=4)
    fig3 = trim_white(fig3, threshold=248, pad=4)

    target_w = max(fig2.width, fig3.width)
    if fig2.width != target_w:
        fig2 = fig2.resize((target_w, round(fig2.height * target_w / fig2.width)), Image.LANCZOS)
    if fig3.width != target_w:
        fig3 = fig3.resize((target_w, round(fig3.height * target_w / fig3.width)), Image.LANCZOS)

    gap = 28
    pad = 18

    canvas_w = target_w + pad * 2
    canvas_h = fig2.height + fig3.height + gap + pad * 2
    canvas = Image.new("RGB", (canvas_w, canvas_h), (255, 255, 255))

    x2 = (canvas_w - fig2.width) // 2
    y2 = pad
    x3 = (canvas_w - fig3.width) // 2
    y3 = y2 + fig2.height + gap

    canvas.paste(fig2, (x2, y2))
    canvas.paste(fig3, (x3, y3))

    canvas.save(out_file)
    tif_file = out_file.rsplit(".", 1)[0] + ".tif"
    canvas.save(tif_file, compression=None)

    print("Saved:")
    print(out_file)
    print(tif_file)

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Stack two graph-summary images vertically into one composite panel.")
    parser.add_argument("--top-image", required=True, help="Path to the upper image.")
    parser.add_argument("--bottom-image", required=True, help="Path to the lower image.")
    parser.add_argument("--output", required=True, help="Output image path.")
    args = parser.parse_args()
    stack_graph_summary_panels(args.top_image, args.bottom_image, args.output)

if __name__ == "__main__":
    main()
