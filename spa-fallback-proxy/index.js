const { Storage } = require('@google-cloud/storage');
const express = require('express');
const app = express();

const BUCKET_NAME = process.env.BUCKET_NAME; // Set this in your environment
const storage = new Storage();

app.get('*', async (req, res) => {
  const filePath = req.path === '/' ? '/index.html' : req.path;
  const file = storage.bucket(BUCKET_NAME).file(filePath.slice(1)); // remove leading '/'

  try {
    await file.exists().then(async ([exists]) => {
      if (exists) {
        file.createReadStream().pipe(res);
      } else {
        // Fallback to index.html for SPA
        storage.bucket(BUCKET_NAME).file('index.html').createReadStream().pipe(res);
      }
    });
  } catch (err) {
    res.status(500).send('Server error');
  }
});

const PORT = process.env.PORT || 8080;
app.listen(PORT, () => {
  console.log(`SPA fallback proxy listening on port ${PORT}`);
});

module.exports = app; 