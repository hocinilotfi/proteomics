
# Style my picture Web Server Gateway Interface

from msproteomics import app
application = app.server
if __name__ == '__main__':
    #app.run(host='127.0.0.1', port=8000, debug=True)
    
    application.run(debug=False)
