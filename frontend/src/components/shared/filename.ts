export function middleEllipsis(filename: string, maxLength: number = 32): string {
  if (filename.length <= maxLength) return filename;
  const extIndex = filename.lastIndexOf('.');
  if (extIndex === -1 || extIndex === 0) return filename.slice(0, maxLength - 3) + '...';
  const ext = filename.slice(extIndex);
  const base = filename.slice(0, maxLength - ext.length - 3);
  return base + '...' + ext;
} 