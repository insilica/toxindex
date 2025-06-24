export async function uploadFileToBackend(file: File, environmentId: string) {
  const formData = new FormData();
  formData.append('file', file);
  formData.append('environment_id', environmentId);
  const res = await fetch('/api/upload-file', {
    method: 'POST',
    credentials: 'include',
    cache: 'no-store',
    body: formData,
  });
  if (!res.ok) throw new Error('Upload failed');
  return await res.json();
} 