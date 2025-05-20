# ===  POS Interface (Backend Only) ===
# Phat Nguyen Cong - Github: https://github.com/paht2005


import pandas as pd

# Simulate Product DB
PRODUCTS = {
    "001": {"name": "Paracetamol", "price": 1.5},
    "002": {"name": "Ibuprofen", "price": 2.0},
    "003": {"name": "ORS Pack", "price": 0.8}
}

# Create invoice
def create_invoice(cart_items: dict):
    lines = []
    total = 0
    for pid, qty in cart_items.items():
        product = PRODUCTS.get(pid)
        if product:
            subtotal = product['price'] * qty
            lines.append({"product": product['name'], "qty": qty, "unit": product['price'], "subtotal": subtotal})
            total += subtotal
    return lines, total
