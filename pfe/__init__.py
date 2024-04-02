from flask import Flask

app = Flask(__name__)


@app.after_request
def add_header(response):
    response.cache_control.private = True
    return response


