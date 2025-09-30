---
title: redeploy docusaurus
sidebar_position: 1
---

cd guide && npm run build
gsutil -m rsync -r build/guide gs://toxindex-react/guide
gsutil -m setmeta -h "Cache-Control:no-store" gs://toxindex-react/**
gcloud compute url-maps invalidate-cdn-cache frontend-url-map  --path "/*"