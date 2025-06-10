from flask import Flask
app = Flask(__name__)

@app.route("/api/me")
def api_me():
    return "ALIVE"

if __name__ == "__main__":
    app.run(port=6513, debug=True) 