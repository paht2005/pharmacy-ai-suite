# === Barcode Scan and Price Tag Detection Tool ===
# Author: Phat Nguyen Cong
# GitHub: https://github.com/paht2005

import cv2
import easyocr
import matplotlib.pyplot as plt
import os
import argparse
import logging
from typing import List, Tuple

logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(message)s')

def load_image(image_path: str):
    if not os.path.exists(image_path):
        logging.error(f"Image path '{image_path}' does not exist.")
        raise FileNotFoundError(f"File not found: {image_path}")
    image = cv2.imread(image_path)
    if image is None:
        logging.error(f"Failed to load image: {image_path}")
        raise ValueError(f"Could not read the image at {image_path}")
    return image

def is_price(text: str) -> bool:
    text = text.strip().replace(',', '')
    return text.startswith('$') or (
        text.count('.') == 1 and text.replace('.', '').isdigit() and 1 < len(text) < 10
    )

def annotate_prices(image, detections: List[Tuple], prices: List[str]):
    for (bbox, text, _) in detections:
        if text in prices:
            top_left = tuple(map(int, bbox[0]))
            bottom_right = tuple(map(int, bbox[2]))
            cv2.rectangle(image, top_left, bottom_right, (0, 255, 0), 2)
            cv2.putText(image, text, (top_left[0], top_left[1] - 10),
                        cv2.FONT_HERSHEY_SIMPLEX, 0.6, (0, 0, 255), 2)

def detect_prices(image_path: str) -> List[str]:
    image = load_image(image_path)
    rgb_image = cv2.cvtColor(image, cv2.COLOR_BGR2RGB)

    reader = easyocr.Reader(['en'], gpu=False)
    detections = reader.readtext(image)

    prices = [text for (_, text, _) in detections if is_price(text)]
    annotate_prices(rgb_image, detections, prices)

    plt.figure(figsize=(12, 8))
    plt.imshow(rgb_image)
    plt.axis('off')
    plt.title("Detected Price Tags")
    plt.tight_layout()
    plt.show()

    return prices

def main():
    parser = argparse.ArgumentParser(description="Detect price tags in an image using OCR.")
    parser.add_argument("image_path", nargs='?', default="test_image.jpg", help="Path to the image file")
    args = parser.parse_args()

    try:
        prices = detect_prices(args.image_path)
        logging.info(f"Detected prices: {prices if prices else 'None found'}")
    except Exception as e:
        logging.error(str(e))

if __name__ == "__main__":
    main()
