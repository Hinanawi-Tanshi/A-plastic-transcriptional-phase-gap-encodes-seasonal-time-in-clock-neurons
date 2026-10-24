function exists = bootstrapRawIsFolder(locations)
exists = false(size(locations));
for k = 1:numel(locations), exists(k) = bootstrapRawStore('isfolder', locations(k)); end
end
