# === Barcode Scan and Price Check ===
# Phat Nguyen Cong 
#  Github: https://github.com/paht2005

import cv2
import easyocr
import matplotlib.pyplot as plt

def detect_prices_from_image(image_path):
    image = cv2.imread(image_path)
    image_rgb = cv2.cvtColor(image, cv2.COLOR_BGR2RGB)
    reader = easyocr.Reader(['en'])
    results = reader.readtext(image)

    prices = []
    for (bbox, text, prob) in results:
        if "$" in text or text.replace(".", "").isdigit():
            prices.append(text)
            top_left = tuple(map(int, bbox[0]))
            bottom_right = tuple(map(int, bbox[2]))
            cv2.rectangle(image_rgb, top_left, bottom_right, (0, 255, 0), 2)
            cv2.putText(image_rgb, text, top_left, cv2.FONT_HERSHEY_SIMPLEX, 0.7, (255, 0, 0), 2)

    plt.imshow(image_rgb)
    plt.axis('off')
    plt.title("🧾Detected Price Tags")
    plt.show()
    return prices

# Example Test
if __name__ == "__main__":
    print(1)
    # Test OCR image 
    # detect_prices_from_image("store_shelf.jpg")
